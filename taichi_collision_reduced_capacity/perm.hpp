// P1920 permutation (the real TaiChi permutation) with a verified inverse.
// Forward matches the NGCC reference bit-for-bit (checked in perm_test).
#pragma once
#include <cstdint>
#include <cstring>

namespace taichi {
typedef uint64_t u64;

static inline u64 ROTL(u64 x, unsigned n) { return n ? ((x << n) | (x >> (64 - n))) : x; }
static inline u64 ROTR(u64 x, unsigned n) { return n ? ((x >> n) | (x << (64 - n))) : x; }

// forward straight-line recurrence tables (from the reference)
static const uint8_t N1[16] = {0, 1, 2, 4, 3, 5, 6, 7, 13, 14, 15, 16, 15, 18, 11, 20};
static const uint8_t N2[16] = {1, 2, 4, 3, 5, 6, 0, 11, 8, 9, 10, 12, 17, 14, 19, 13};
static const uint8_t MM[16] = {8, 0, 0, 0, 0, 0, 23, 62, 35, 14, 48, 1, 57, 63, 58, 22};
static const unsigned OROT[6] = {0, 15, 48, 63, 6, 33};  // output rotations of the MRM

struct Tables {
  u64 RC[120];
  u64 Hinv;                 // ring-inverse mask of H(w)=w^ROTL(w,23)^ROTL(w,31)
  uint32_t sinv_anf[5];     // ANF monomial masks of the inverse S-box (over b0..b4)
  int k2fwd[30];            // forward lane permutation: next[k2fwd[n]] = cur[n]

  Tables() {
    build_RC();
    build_k2();
    build_sbox_inverse();
    build_Hinv();
  }
  void build_RC() {
    static const int p[7] = {0, 2, 3, 6, 13, 28, 59};
    uint8_t s[7] = {1, 0, 0, 0, 0, 0, 0};
    for (int q = 0; q < 120; q++) {
      u64 C = 0;
      for (int i = 0; i < 7; i++) if (s[i]) C ^= (1ULL << p[i]);
      RC[q] = C;
      uint8_t n[7];
      n[0] = s[6]; n[1] = s[0] ^ s[6]; n[2] = s[1]; n[3] = s[2]; n[4] = s[3]; n[5] = s[4]; n[6] = s[5];
      for (int i = 0; i < 7; i++) s[i] = n[i];
    }
  }
  void build_k2() {
    for (int n = 0; n < 30; n++) k2fwd[n] = (10 * (n + 1) % 31) - 1;
  }
  // forward 5-bit S-box on (a0..a4), bit j = a_j
  static uint32_t sbox_fwd_scalar(uint32_t v) {
    uint32_t a0 = v & 1, a1 = (v >> 1) & 1, a2 = (v >> 2) & 1, a3 = (v >> 3) & 1, a4 = (v >> 4) & 1;
    uint32_t b0 = a1 ^ (a0 & a1) ^ a2 ^ (a1 & a2) ^ a3 ^ (a3 & a4);
    uint32_t b1 = 1u ^ a1 ^ (a0 & a3) ^ (a1 & a3) ^ a4 ^ (a2 & a4);
    uint32_t b2 = (a1 & a2) ^ a3 ^ (a2 & a3) ^ a4 ^ (a0 & a4);
    uint32_t b3 = (a0 & a2) ^ (a1 & a3) ^ (a2 & a3) ^ a4;
    uint32_t b4 = a0 ^ (a2 & a3) ^ (a1 & a4);
    return b0 | (b1 << 1) | (b2 << 2) | (b3 << 3) | (b4 << 4);
  }
  bool sbox_is_bijection() const {
    int seen[32] = {0};
    for (uint32_t v = 0; v < 32; v++) seen[sbox_fwd_scalar(v)]++;
    for (int i = 0; i < 32; i++) if (seen[i] != 1) return false;
    return true;
  }
  void build_sbox_inverse() {
    uint32_t sinv[32];
    for (uint32_t v = 0; v < 32; v++) sinv[sbox_fwd_scalar(v)] = v;
    // For each output bit k, truth table over input index w (= b0..b4): bit k of sinv[w].
    for (int k = 0; k < 5; k++) {
      uint32_t f[32];
      for (int w = 0; w < 32; w++) f[w] = (sinv[w] >> k) & 1;
      // Mobius transform -> ANF coefficients
      for (int i = 0; i < 5; i++)
        for (int w = 0; w < 32; w++)
          if (w & (1 << i)) f[w] ^= f[w ^ (1 << i)];
      uint32_t mask = 0;
      for (int m = 0; m < 32; m++) if (f[m]) mask |= (1u << m);
      sinv_anf[k] = mask;
    }
  }
  void build_Hinv() {
    // H is multiplication by (1 + x^23 + x^31) in GF(2)[x]/(x^64+1). Invert the 64x64
    // circulant by Gaussian elimination, then Hinv = inverse applied to the unit 1.
    u64 M[64];
    for (int i = 0; i < 64; i++) {
      u64 e = 1ULL << i;
      M[i] = e ^ ROTL(e, 23) ^ ROTL(e, 31);   // column i of H as a 64-bit vector
    }
    // Treat M[i] as the image of basis vector e_i; solve M * h = e0 by reducing columns.
    // Build row-wise matrix A where A[r] has bit c = (M[c]>>r)&1; invert A, apply to e0.
    u64 A[64], B[64];
    for (int r = 0; r < 64; r++) {
      u64 row = 0;
      for (int c = 0; c < 64; c++) if ((M[c] >> r) & 1) row |= (1ULL << c);
      A[r] = row;
      B[r] = 1ULL << r;
    }
    for (int col = 0; col < 64; col++) {
      int piv = -1;
      for (int r = col; r < 64; r++) if ((A[r] >> col) & 1) { piv = r; break; }
      u64 ta = A[col]; A[col] = A[piv]; A[piv] = ta;
      u64 tb = B[col]; B[col] = B[piv]; B[piv] = tb;
      for (int r = 0; r < 64; r++)
        if (r != col && ((A[r] >> col) & 1)) { A[r] ^= A[col]; B[r] ^= B[col]; }
    }
    // Now B is A^{-1}. h = A^{-1} * e0 : bit r of h = (B[r]>>0)&1.
    u64 h = 0;
    for (int r = 0; r < 64; r++) if (B[r] & 1ULL) h |= (1ULL << r);
    Hinv = h;
  }
  inline u64 apply_Hinv(u64 v) const {
    u64 r = 0, h = Hinv;
    while (h) { int j = __builtin_ctzll(h); r ^= ROTL(v, j); h &= h - 1; }
    return r;
  }
};

extern const Tables T;

// ---- forward P1920 (matches reference) ----
inline void P1920_fwd(u64 A[30], int b) {
  u64 tmp[30];
  u64 *cur = A, *next = tmp;
  for (int t = 0; t < 12; t++) {
    for (int j = 0; j < 5; j++) {
      u64 y[22];
      for (int i = 0; i < 6; i++) y[i] = cur[5 * i + j];
      for (int i = 6; i < 22; i++) y[i] = y[N1[i - 6]] ^ ROTL(y[N2[i - 6]], MM[i - 6]);
      for (int i = 0; i < 6; i++) cur[5 * i + j] = ROTL(y[16 + i], OROT[i]);
    }
    for (int i = 0; i < 6; i++) {
      u64 a0 = cur[5 * i + 0], a1 = cur[5 * i + 1], a2 = cur[5 * i + 2], a3 = cur[5 * i + 3], a4 = cur[5 * i + 4];
      cur[5 * i + 0] = a1 ^ (a0 & a1) ^ a2 ^ (a1 & a2) ^ a3 ^ (a3 & a4);
      cur[5 * i + 1] = ~0ULL ^ a1 ^ (a0 & a3) ^ (a1 & a3) ^ a4 ^ (a2 & a4);
      cur[5 * i + 2] = (a1 & a2) ^ a3 ^ (a2 & a3) ^ a4 ^ (a0 & a4);
      cur[5 * i + 3] = (a0 & a2) ^ (a1 & a3) ^ (a2 & a3) ^ a4;
      cur[5 * i + 4] = a0 ^ (a2 & a3) ^ (a1 & a4);
    }
    for (int n = 0; n < 30; n++) next[T.k2fwd[n]] = cur[n];
    for (int j = 0; j < 5; j++) next[j] ^= T.RC[(b - 1) * 60 + 5 * t + j];
    u64 *sw = cur; cur = next; next = sw;
  }
  if (cur != A) std::memcpy(A, tmp, 30 * sizeof(u64));
}

// ---- inverse layers ----
// H(w) = w ^ ROTL(w,23) ^ ROTL(w,31) is multiplication by f = 1 + x^23 + x^31 in
// GF(2)[x]/(x^64+1). Frobenius gives f^8 = 1 + x^184 + x^248 = 1 + 2x^56 = 1, so
// f^{-1} = f^7 = f * f^2 * f^4 with f^2 = 1 + x^46 + x^62 and f^4 = 1 + x^28 + x^60.
static inline u64 Hinv_fast(u64 v) {
  v = v ^ ROTL(v, 23) ^ ROTL(v, 31);
  v = v ^ ROTL(v, 46) ^ ROTL(v, 62);
  v = v ^ ROTL(v, 28) ^ ROTL(v, 60);
  return v;
}

// Inverse S-box, explicit degree-3 ANF (computed by Mobius transform of the inverted
// 5-bit table; cross-checked against the table version in the tests).
static inline void sbox_inv5(u64& r0, u64& r1, u64& r2, u64& r3, u64& r4) {
  const u64 b0 = r0, b1 = r1, b2 = r2, b3 = r3, b4 = r4;
  const u64 b01 = b0 & b1, b02 = b0 & b2, b03 = b0 & b3, b04 = b0 & b4, b12 = b1 & b2;
  const u64 b13 = b1 & b3, b14 = b1 & b4, b23 = b2 & b3, b24 = b2 & b4, b34 = b3 & b4;
  const u64 b012 = b01 & b2, b013 = b01 & b3, b014 = b01 & b4, b023 = b02 & b3, b024 = b02 & b4;
  const u64 b034 = b03 & b4, b123 = b12 & b3, b124 = b12 & b4, b134 = b13 & b4, b234 = b23 & b4;
  r0 = ~(b0 ^ b1 ^ b01 ^ b12 ^ b012 ^ b3 ^ b013 ^ b14 ^ b024 ^ b34 ^ b134);
  r1 = b0 ^ b01 ^ b2 ^ b02 ^ b03 ^ b13 ^ b013 ^ b23 ^ b4 ^ b04 ^ b14 ^ b014 ^ b24 ^ b124 ^ b34;
  r2 = ~(b0 ^ b1 ^ b02 ^ b3 ^ b13 ^ b023 ^ b4 ^ b04 ^ b14 ^ b014 ^ b24 ^ b34 ^ b034);
  r3 = ~(b0 ^ b1 ^ b01 ^ b2 ^ b02 ^ b012 ^ b03 ^ b4 ^ b14 ^ b134 ^ b234);
  r4 = b2 ^ b02 ^ b3 ^ b23 ^ b123 ^ b04 ^ b024 ^ b124 ^ b034 ^ b134 ^ b234;
}

// Table-driven inverse S-box (slow; test oracle only).
inline void sbox_inv_row_table(u64 r[5]) {
  // r currently holds (b0..b4); produce (a0..a4) via the inverse ANF over b.
  u64 P[32];
  P[0] = ~0ULL;
  for (int m = 1; m < 32; m++) { int lb = __builtin_ctz(m); P[m] = P[m ^ (1 << lb)] & r[lb]; }
  u64 out[5];
  for (int k = 0; k < 5; k++) {
    u64 acc = 0;
    uint32_t mask = T.sinv_anf[k];
    while (mask) { int m = __builtin_ctz(mask); acc ^= P[m]; mask &= mask - 1; }
    out[k] = acc;
  }
  for (int k = 0; k < 5; k++) r[k] = out[k];
}

inline void mrm_inv_col(u64 o[6]) {
  u64 y16 = o[0], y17 = ROTR(o[1], 15), y18 = ROTR(o[2], 48), y19 = ROTR(o[3], 63),
      y20 = ROTR(o[4], 6), y21 = ROTR(o[5], 33);
  u64 y13 = ROTR(y21 ^ y20, 22);
  u64 y11 = y20 ^ ROTL(y19, 58);
  u64 y14 = ROTR(y19 ^ y18, 63);
  u64 y15 = y18 ^ ROTL(y17, 57);
  u64 y12 = ROTR(y17 ^ y16, 1);
  u64 y10 = ROTR(y16 ^ y15, 48);
  u64 y9 = ROTR(y15 ^ y14, 14);
  u64 y8 = ROTR(y14 ^ y13, 35);
  u64 y7 = y13 ^ ROTL(y11, 62);
  u64 K = y7 ^ y8 ^ y9 ^ y10 ^ y11;
  u64 rhs = y12 ^ K ^ ROTL(K, 23);
  u64 y1 = Hinv_fast(rhs);
  u64 y2 = y7 ^ y1;
  u64 y4 = y8 ^ y2;
  u64 y3 = y9 ^ y4;
  u64 y5 = y10 ^ y3;
  u64 y6 = K ^ y1;
  u64 y0 = y6 ^ ROTL(y1, 8);
  o[0] = y0; o[1] = y1; o[2] = y2; o[3] = y3; o[4] = y4; o[5] = y5;
}

inline void P1920_inv(u64 A[30], int b) {
  for (int t = 11; t >= 0; t--) {
    // undo RC
    for (int j = 0; j < 5; j++) A[j] ^= T.RC[(b - 1) * 60 + 5 * t + j];
    // undo K2 lane permutation: pre[n] = post[k2fwd[n]]
    u64 pre[30];
    for (int n = 0; n < 30; n++) pre[n] = A[T.k2fwd[n]];
    // undo S-box per row
    for (int i = 0; i < 6; i++)
      sbox_inv5(pre[5 * i + 0], pre[5 * i + 1], pre[5 * i + 2], pre[5 * i + 3], pre[5 * i + 4]);
    // undo MRM per column
    for (int j = 0; j < 5; j++) {
      u64 o[6];
      for (int i = 0; i < 6; i++) o[i] = pre[5 * i + j];
      mrm_inv_col(o);
      for (int i = 0; i < 6; i++) pre[5 * i + j] = o[i];
    }
    std::memcpy(A, pre, 30 * sizeof(u64));
  }
}

// ---- batched forward/inverse: BS independent states, lane-minor layout st[30][BS] ----
static const int BS = 64;

inline void P1920_fwd_batch(u64 st[30][BS], int b) {
  alignas(64) u64 tmp[30][BS];
  u64 (*cur)[BS] = st, (*next)[BS] = tmp;
  for (int t = 0; t < 12; t++) {
    for (int j = 0; j < 5; j++) {
      alignas(64) u64 y[22][BS];
      for (int i = 0; i < 6; i++)
        for (int l = 0; l < BS; l++) y[i][l] = cur[5 * i + j][l];
      for (int i = 6; i < 22; i++) {
        const int a = N1[i - 6], c = N2[i - 6], mm = MM[i - 6];
        if (mm) for (int l = 0; l < BS; l++) y[i][l] = y[a][l] ^ ROTL(y[c][l], mm);
        else    for (int l = 0; l < BS; l++) y[i][l] = y[a][l] ^ y[c][l];
      }
      for (int i = 0; i < 6; i++) {
        const unsigned rr = OROT[i];
        if (rr) for (int l = 0; l < BS; l++) cur[5 * i + j][l] = ROTL(y[16 + i][l], rr);
        else    for (int l = 0; l < BS; l++) cur[5 * i + j][l] = y[16 + i][l];
      }
    }
    for (int i = 0; i < 6; i++) {
      u64 *a0 = cur[5 * i + 0], *a1 = cur[5 * i + 1], *a2 = cur[5 * i + 2], *a3 = cur[5 * i + 3], *a4 = cur[5 * i + 4];
      for (int l = 0; l < BS; l++) {
        u64 x0 = a0[l], x1 = a1[l], x2 = a2[l], x3 = a3[l], x4 = a4[l];
        a0[l] = x1 ^ (x0 & x1) ^ x2 ^ (x1 & x2) ^ x3 ^ (x3 & x4);
        a1[l] = ~0ULL ^ x1 ^ (x0 & x3) ^ (x1 & x3) ^ x4 ^ (x2 & x4);
        a2[l] = (x1 & x2) ^ x3 ^ (x2 & x3) ^ x4 ^ (x0 & x4);
        a3[l] = (x0 & x2) ^ (x1 & x3) ^ (x2 & x3) ^ x4;
        a4[l] = x0 ^ (x2 & x3) ^ (x1 & x4);
      }
    }
    for (int n = 0; n < 30; n++) std::memcpy(next[T.k2fwd[n]], cur[n], BS * sizeof(u64));
    for (int j = 0; j < 5; j++) { const u64 rc = T.RC[(b - 1) * 60 + 5 * t + j]; for (int l = 0; l < BS; l++) next[j][l] ^= rc; }
    u64 (*sw)[BS] = cur; cur = next; next = sw;
  }
  if (cur != st) std::memcpy(st, tmp, 30 * BS * sizeof(u64));
}

inline void P1920_inv_batch(u64 st[30][BS], int b) {
  for (int t = 11; t >= 0; t--) {
    for (int j = 0; j < 5; j++) { const u64 rc = T.RC[(b - 1) * 60 + 5 * t + j]; for (int l = 0; l < BS; l++) st[j][l] ^= rc; }
    alignas(64) u64 pre[30][BS];
    for (int n = 0; n < 30; n++) std::memcpy(pre[n], st[T.k2fwd[n]], BS * sizeof(u64));
    for (int i = 0; i < 6; i++) {
      u64 *r0 = pre[5 * i + 0], *r1 = pre[5 * i + 1], *r2 = pre[5 * i + 2], *r3 = pre[5 * i + 3], *r4 = pre[5 * i + 4];
      for (int l = 0; l < BS; l++) sbox_inv5(r0[l], r1[l], r2[l], r3[l], r4[l]);
    }
    for (int j = 0; j < 5; j++) {
      u64 *c0 = pre[0 * 5 + j], *c1 = pre[1 * 5 + j], *c2 = pre[2 * 5 + j], *c3 = pre[3 * 5 + j], *c4 = pre[4 * 5 + j], *c5 = pre[5 * 5 + j];
      for (int l = 0; l < BS; l++) {
        u64 y16 = c0[l], y17 = ROTR(c1[l], 15), y18 = ROTR(c2[l], 48), y19 = ROTR(c3[l], 63), y20 = ROTR(c4[l], 6), y21 = ROTR(c5[l], 33);
        u64 y13 = ROTR(y21 ^ y20, 22), y11 = y20 ^ ROTL(y19, 58), y14 = ROTR(y19 ^ y18, 63), y15 = y18 ^ ROTL(y17, 57), y12 = ROTR(y17 ^ y16, 1);
        u64 y10 = ROTR(y16 ^ y15, 48), y9 = ROTR(y15 ^ y14, 14), y8 = ROTR(y14 ^ y13, 35), y7 = y13 ^ ROTL(y11, 62);
        u64 K = y7 ^ y8 ^ y9 ^ y10 ^ y11;
        u64 y1 = Hinv_fast(y12 ^ K ^ ROTL(K, 23));
        u64 y2 = y7 ^ y1, y4 = y8 ^ y2, y3 = y9 ^ y4, y5 = y10 ^ y3, y6 = K ^ y1, y0 = y6 ^ ROTL(y1, 8);
        c0[l] = y0; c1[l] = y1; c2[l] = y2; c3[l] = y3; c4[l] = y4; c5[l] = y5;
      }
    }
    std::memcpy(st, pre, 30 * BS * sizeof(u64));
  }
}

}  // namespace taichi
