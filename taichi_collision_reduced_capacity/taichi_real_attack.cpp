// taichi_real_attack.cpp
//
// Structural 2^{3c/4} collision attack on the TaiChi IFS mode, run against the REAL TaiChi
// permutation P1920 (P1 = round-constant bank 1, P2 = bank 2) with the capacity reduced to
// c bits and the rate set to 1920 - c bits, as in any sponge.
//
// Reduced layout (c a multiple of 8, 8 <= c <= 64):
//   * rate     = the first (1920-c)/8 bytes of the 240-byte state (memcpy of RATE_BYTES as in
//                the reference), i.e. words 0..28 and the low 64-c bits of word 29;
//   * capacity = the last c/8 bytes = the high c bits of word 29 (L and R are c-bit values);
//   * domain   = tau_Abs = 1 / tau_Fin = 2, XORed into the low bits of the capacity window
//                (the reference XORs it into the low byte of word 29, which is capacity there);
//   * 10*1 padding exactly as in the reference, with RATE_BITS = 1920 - c;
//   * digest   = first c/8 bytes of the final rate S (the collision is on all of S anyway).
// With c = 512 and the capacity in words 22..29 this is the original TaiChi-512; the selftest
// reproduces the reference digests in that configuration.
//
// Build: g++ -O3 -march=native -std=c++17 -pthread taichi_real_attack.cpp -o taichi_real
//        (x86 with AVX-512: add -mprefer-vector-width=512; Apple silicon: -mcpu=native)
//
// Modes:
//   selftest                      permutation, inverse, real TaiChi-512 digests, padding, engines
//   bench                         permutation throughput (1 and T threads)
//   single --c 24                 interleaved run (c <= 32), stop at the first D-collision, report
//   scale  --cs 8,16,24,32        first-collision q over many independent A*, m0, z choices
//   batch  --c 40                 large c: memory-planned join, stops at the first D-collision
// Options: --seed S --threads T --trials N --mem-gb G --lambda L --q Q --qf Q --qi Q --no-stop

#include "perm.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <sys/mman.h>

namespace taichi { const Tables T; }
using namespace taichi;
typedef uint32_t u32;
typedef unsigned __int128 u128;

static const u64 TAU_ABS = 1;
static const u64 TAU_FIN = 2;

static double now_sec() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static inline u64 mix64(u64 z) {
  z ^= z >> 30; z *= 0xbf58476d1ce4e5b9ULL;
  z ^= z >> 27; z *= 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
struct SplitMix64 {
  u64 s;
  explicit SplitMix64(u64 seed) : s(seed) {}
  u64 next() { s += 0x9e3779b97f4a7c15ULL; return mix64(s); }
};
static std::string hexw(u64 v, int bits) {
  char buf[40];
  snprintf(buf, sizeof buf, "0x%0*" PRIx64, (bits + 3) / 4, v);
  return buf;
}
static std::string hexbytes(const uint8_t* p, size_t n) {
  std::string s;
  s.reserve(2 * n);
  static const char* hx = "0123456789abcdef";
  for (size_t k = 0; k < n; k++) { s.push_back(hx[p[k] >> 4]); s.push_back(hx[p[k] & 15]); }
  return s;
}
static std::string hexwords(const u64* w, int n) {
  std::string s;
  char buf[24];
  for (int k = 0; k < n; k++) { snprintf(buf, sizeof buf, "%s%016" PRIx64, k ? " " : "", w[k]); s += buf; }
  return s;
}
#define HX(v, b) hexw((v), (b)).c_str()

// =========================================================================================
// Reference padding (exact port of TaiChi_Padding, RATE_BITS as a parameter).
// =========================================================================================
static u64 ref_padding(std::vector<uint8_t>& out, const uint8_t* msg, u64 msg_len_bits, u64 RATE_BITS) {
  const u64 RATE_BYTES = RATE_BITS / 8;
  const u64 padded_len = (msg_len_bits + 2 + RATE_BITS - 1) / RATE_BITS * RATE_BYTES;
  out.assign(padded_len + 8, 0);
  const u64 bytes = msg_len_bits / 8;
  if (bytes) std::memcpy(out.data(), msg, bytes);
  const u64 remain = msg_len_bits % 8;
  if (remain) out[bytes] = (uint8_t)(msg[bytes] & (0xff << (8 - remain)));
  u64 total = msg_len_bits;
  out[total / 8] |= (uint8_t)(1 << (7 - total % 8));
  total++;
  while ((total + 1) % RATE_BITS != 0) total++;
  out[total / 8] |= (uint8_t)(1 << (7 - total % 8));
  total++;
  return (total + 7) / 8;
}

// =========================================================================================
// Real TaiChi-512 (rate 1408 = words 0..21, capacity = words 22..29, domain in word 29).
// Used only to show that P1920 and the mode wiring reproduce the reference digests.
// =========================================================================================
static void real_taichi512(const uint8_t* msg, u64 bits, uint8_t dig[64]) {
  const int RW = 22, CW = 8;
  const u64 RATE_BITS = 1408, RB = 176;
  std::vector<uint8_t> buf;
  const u64 plen = ref_padding(buf, msg, bits, RATE_BITS);
  u64 S[22] = {0}, L[8] = {0}, R[8] = {0};
  const u64 nb = plen / RB;
  for (u64 i = 0; i < nb; i++) {
    u64 M[22] = {0};
    std::memcpy(M, buf.data() + i * RB, RB);
    u64 st[30] = {0};
    for (int k = 0; k < RW; k++) st[k] = S[k] ^ M[k];
    for (int k = 0; k < CW; k++) st[RW + k] = L[k];
    st[29] ^= (i + 1 == nb) ? TAU_FIN : TAU_ABS;
    P1920_fwd(st, 1);
    u64 B[8];
    for (int k = 0; k < CW; k++) B[k] = st[RW + k];
    for (int k = 0; k < CW; k++) st[RW + k] = R[k] ^ B[k];
    P1920_fwd(st, 2);
    for (int k = 0; k < RW; k++) S[k] = st[k];
    for (int k = 0; k < CW; k++) { L[k] = st[RW + k]; R[k] = B[k] ^ st[RW + k]; }
  }
  std::memcpy(dig, S, 64);
}

// =========================================================================================
// Reduced layout: capacity = high c bits of word 29.
// =========================================================================================
struct Layout {
  int c = 40;
  u64 R = 0, RB = 0;          // rate bits / bytes
  int capshift = 0;            // 64 - c
  u64 lowmask = 0, cmask = 0;  // rate bits of word 29 / c-bit mask
  int pad_word = 0, bitA = 0, bitB = 0;
  u64 padmask = 0;
  void init(int c_) {
    c = c_;
    R = 1920 - (u64)c;
    RB = R / 8;
    capshift = 64 - c;
    lowmask = capshift ? ((1ULL << capshift) - 1) : 0;
    cmask = (c == 64) ? ~0ULL : ((1ULL << c) - 1);
    // The last two bits of a full final block (positions R-2, R-1) are the "11" of 10*1 when
    // the message has 2R-2 bits. They sit in byte RB-1, bit masks 0x02 / 0x01 (MSB-first).
    const u64 lb = RB - 1;
    pad_word = (int)(lb / 8);
    bitA = (int)((lb % 8) * 8 + 0);
    bitB = (int)((lb % 8) * 8 + 1);
    padmask = (1ULL << bitA) | (1ULL << bitB);
  }
  inline u64 cap_get(u64 w29) const { return capshift ? (w29 >> capshift) : w29; }
  inline u64 cap_set(u64 w29, u64 v) const {
    return capshift ? ((w29 & lowmask) | ((v & cmask) << capshift)) : v;
  }
};

// Reduced reference hash on raw bytes: the reference code with RATE_BITS = 1920 - c.
struct RTrace {
  std::vector<std::array<u64, 30>> blocks;
  u64 S[30] = {0};
  u64 L = 0, Rr = 0;
  u64 p1in[30] = {0}, p2in[30] = {0};   // inputs of the last P1 / P2 call
  std::vector<uint8_t> digest;
};
static RTrace reduced_hash(const Layout& Ly, const uint8_t* msg, u64 bits) {
  RTrace tr;
  std::vector<uint8_t> buf;
  const u64 plen = ref_padding(buf, msg, bits, Ly.R);
  const u64 nb = plen / Ly.RB;
  u64 S[30] = {0}, Lc = 0, Rc = 0;
  for (u64 i = 0; i < nb; i++) {
    std::array<u64, 30> M{};
    std::memcpy(M.data(), buf.data() + i * Ly.RB, Ly.RB);
    tr.blocks.push_back(M);
    const u64 tau = (i + 1 == nb) ? TAU_FIN : TAU_ABS;
    u64 st[30];
    for (int k = 0; k < 30; k++) st[k] = S[k] ^ M[k];
    st[29] = Ly.cap_set(st[29], Lc ^ tau);
    std::memcpy(tr.p1in, st, sizeof st);
    P1920_fwd(st, 1);
    const u64 B = Ly.cap_get(st[29]);
    st[29] = Ly.cap_set(st[29], Rc ^ B);
    std::memcpy(tr.p2in, st, sizeof st);
    P1920_fwd(st, 2);
    const u64 C = Ly.cap_get(st[29]);
    std::memcpy(S, st, sizeof st);
    S[29] &= Ly.lowmask;
    Lc = C;
    Rc = (B ^ C) & Ly.cmask;
  }
  std::memcpy(tr.S, S, sizeof S);
  tr.L = Lc;
  tr.Rr = Rc;
  tr.digest.assign((const uint8_t*)S, (const uint8_t*)S + Ly.c / 8);
  return tr;
}

// =========================================================================================
// Attack instance: fixed real permutation; A*, the m0 family and the z family vary per seed.
// =========================================================================================
struct Instance {
  Layout Ly;
  u64 seed = 0;
  u64 Astar[30], m0base[30];
  u64 zbase = 0;
  void init(int c, u64 seed_) {
    Ly.init(c);
    seed = seed_;
    SplitMix64 g(seed);
    for (int k = 0; k < 30; k++) Astar[k] = g.next();
    Astar[29] &= Ly.lowmask;
    for (int k = 0; k < 30; k++) m0base[k] = g.next();
    m0base[29] &= Ly.lowmask;
    zbase = g.next() & Ly.cmask;
  }
  inline u64 z_of(u64 j) const { return (zbase ^ j) & Ly.cmask; }
  inline void m0_of(u64 i, u64 out[30]) const { std::memcpy(out, m0base, sizeof m0base); out[0] ^= i; }
  // First-block root: P1(m0 || tau_Abs) = (a, x); P2(a || x) = (s, y). Two queries.
  void fwd_root(u64 i, u64 s[30], u64& x, u64& y) const {
    u64 st[30];
    m0_of(i, st);
    st[29] = Ly.cap_set(st[29], TAU_ABS);
    P1920_fwd(st, 1);
    x = Ly.cap_get(st[29]);
    P1920_fwd(st, 2);          // R0 = 0, so the P2 input is exactly the P1 output (a || x)
    y = Ly.cap_get(st[29]);
    std::memcpy(s, st, sizeof st);
    s[29] &= Ly.lowmask;
  }
  // Inverse query: P1^{-1}(A* || z) = (u, l). One query.
  void inv_query(u64 j, u64 u[30], u64& l) const {
    u64 st[30];
    std::memcpy(st, Astar, sizeof st);
    st[29] = Ly.cap_set(st[29], z_of(j));
    P1920_inv(st, 1);
    l = Ly.cap_get(st[29]);
    std::memcpy(u, st, sizeof st);
    u[29] &= Ly.lowmask;
  }
};

// Batched roots for arbitrary indices: y (and optionally x and the padding word of s).
static void roots_idx(const Instance& I, const u64* idx, int m, u64* y, u64* x, u64* sw) {
  alignas(64) u64 st[30][BS];
  for (int w = 0; w < 30; w++)
    for (int l = 0; l < BS; l++) st[w][l] = I.m0base[w];
  for (int l = 0; l < m; l++) st[0][l] ^= idx[l];
  for (int l = 0; l < BS; l++) st[29][l] = I.Ly.cap_set(st[29][l], TAU_ABS);
  P1920_fwd_batch(st, 1);
  if (x) for (int l = 0; l < m; l++) x[l] = I.Ly.cap_get(st[29][l]);
  P1920_fwd_batch(st, 2);
  for (int l = 0; l < m; l++) y[l] = I.Ly.cap_get(st[29][l]);
  if (sw) for (int l = 0; l < m; l++) sw[l] = st[I.Ly.pad_word][l];
}
static void roots_range(const Instance& I, u64 i0, int m, u64* y, u64* x, u64* sw) {
  u64 idx[BS];
  for (int l = 0; l < m; l++) idx[l] = i0 + l;
  roots_idx(I, idx, m, y, x, sw);
}
// Batched inverse queries j0..j0+m-1: key = l ^ tau_Fin and the padding word of u.
static void queries_range(const Instance& I, u64 j0, int m, u64* key, u64* uw) {
  alignas(64) u64 st[30][BS];
  for (int w = 0; w < 30; w++)
    for (int l = 0; l < BS; l++) st[w][l] = I.Astar[w];
  for (int l = 0; l < BS; l++) st[29][l] = I.Ly.cap_set(st[29][l], I.z_of(j0 + l));
  P1920_inv_batch(st, 1);
  for (int l = 0; l < m; l++) {
    key[l] = I.Ly.cap_get(st[29][l]) ^ TAU_FIN;
    if (uw) uw[l] = st[I.Ly.pad_word][l];
  }
}

// =========================================================================================
// Join bookkeeping, explanation and verification.
// =========================================================================================
struct Surv { u64 D, i, j; };

struct JoinInfo {
  u64 i = 0, j = 0;
  u64 m0[30], a[30], s[30], u[30], m1[30];
  u64 x = 0, y = 0, z = 0, l = 0, D = 0;
  bool key_ok = false, pad_ok = false;
};
static JoinInfo explain_join(const Instance& I, u64 i, u64 j) {
  const Layout& Ly = I.Ly;
  JoinInfo J;
  J.i = i; J.j = j;
  I.m0_of(i, J.m0);
  u64 st[30];
  std::memcpy(st, J.m0, sizeof st);
  st[29] = Ly.cap_set(st[29], TAU_ABS);
  P1920_fwd(st, 1);
  J.x = Ly.cap_get(st[29]);
  std::memcpy(J.a, st, sizeof st);
  J.a[29] &= Ly.lowmask;
  P1920_fwd(st, 2);
  J.y = Ly.cap_get(st[29]);
  std::memcpy(J.s, st, sizeof st);
  J.s[29] &= Ly.lowmask;
  J.z = I.z_of(j);
  std::memcpy(st, I.Astar, sizeof st);
  st[29] = Ly.cap_set(st[29], J.z);
  P1920_inv(st, 1);
  J.l = Ly.cap_get(st[29]);
  std::memcpy(J.u, st, sizeof st);
  J.u[29] &= Ly.lowmask;
  J.key_ok = ((J.l ^ TAU_FIN) == J.y);
  for (int k = 0; k < 30; k++) J.m1[k] = J.s[k] ^ J.u[k];
  J.pad_ok = ((J.m1[Ly.pad_word] & Ly.padmask) == Ly.padmask);
  J.D = (J.x ^ J.y ^ J.z) & Ly.cmask;
  return J;
}

// Raw message = block m0 (R bits) followed by the first R-2 bits of block m1 (2R-2 bits).
static std::vector<uint8_t> raw_message(const Layout& Ly, const u64 m0[30], const u64 m1[30], u64& bits) {
  std::vector<uint8_t> out(2 * Ly.RB);
  std::memcpy(out.data(), m0, Ly.RB);
  std::memcpy(out.data() + Ly.RB, m1, Ly.RB);
  out[2 * Ly.RB - 1] &= 0xFC;     // the last two bits are padding, not message
  bits = 2 * Ly.R - 2;
  return out;
}

struct HashCheck {
  RTrace h;
  std::vector<uint8_t> msg;
  u64 bits = 0;
  bool blocks_ok = false, p1_ok = false, p2_ok = false;
};
static HashCheck hash_check(const Instance& I, const JoinInfo& J) {
  const Layout& Ly = I.Ly;
  HashCheck C;
  C.msg = raw_message(Ly, J.m0, J.m1, C.bits);
  C.h = reduced_hash(Ly, C.msg.data(), C.bits);
  C.blocks_ok = C.h.blocks.size() == 2 && std::memcmp(C.h.blocks[0].data(), J.m0, 240) == 0 &&
                std::memcmp(C.h.blocks[1].data(), J.m1, 240) == 0;
  u64 e1[30], e2[30];
  std::memcpy(e1, J.u, sizeof e1);
  e1[29] = Ly.cap_set(e1[29], J.l);
  std::memcpy(e2, I.Astar, sizeof e2);
  e2[29] = Ly.cap_set(e2[29], J.D);
  C.p1_ok = std::memcmp(C.h.p1in, e1, sizeof e1) == 0;
  C.p2_ok = std::memcmp(C.h.p2in, e2, sizeof e2) == 0;
  return C;
}

static void print_join(const Instance& I, int idx, const JoinInfo& J) {
  const int c = I.Ly.c;
  printf("  candidate %d: forward root i = %" PRIu64 ", inverse query j = %" PRIu64 "\n", idx, J.i, J.j);
  printf("    P1(m0 || tau_Abs) capacity x = %s\n", HX(J.x, c));
  printf("    P2(a || x)        capacity y = %s   (state after block 1 = (s, y, x^y), x^y = %s)\n", HX(J.y, c),
         HX(J.x ^ J.y, c));
  printf("    z = %s,  P1^-1(A* || z) capacity l = %s,  l ^ tau_Fin == y: %s\n", HX(J.z, c), HX(J.l, c),
         J.key_ok ? "yes" : "NO");
  printf("    m1 = s ^ u, padding bits (word %d, bits %d,%d) = %d%d\n", I.Ly.pad_word, I.Ly.bitB, I.Ly.bitA,
         (int)((J.m1[I.Ly.pad_word] >> I.Ly.bitB) & 1), (int)((J.m1[I.Ly.pad_word] >> I.Ly.bitA) & 1));
  printf("    x = %s, y = %s, z = %s, D = x^y^z = %s\n", HX(J.x, c), HX(J.y, c), HX(J.z, c), HX(J.D, c));
}

static bool report_collision(const Instance& I, const JoinInfo& A, const JoinInfo& B, const char* outprefix) {
  const Layout& Ly = I.Ly;
  const int c = Ly.c;
  printf("\nCollision candidates (equal D):\n");
  print_join(I, 1, A);
  print_join(I, 2, B);
  printf("\nCommon final P2 input predicted by the attack: A* || D\n");
  printf("  A* (rate, words 0..29, word 29 low %d bits) = %s\n", Ly.capshift, hexwords(I.Astar, 30).c_str());
  printf("  D  (capacity, high %d bits of word 29)      = %s\n", c, HX(A.D, c));
  HashCheck H[2] = {hash_check(I, A), hash_check(I, B)};
  for (int t = 0; t < 2; t++) {
    const HashCheck& h = H[t];
    printf("\nMessage %s: %" PRIu64 " bits (%zu bytes, last byte holds 6 message bits)\n", t ? "M'" : "M", h.bits,
           h.msg.size());
    printf("  raw (hex, MSB-first bytes) = %s\n", hexbytes(h.msg.data(), h.msg.size()).c_str());
    printf("  padded block 1 (rate bytes) = %s\n", hexbytes((const uint8_t*)h.h.blocks[0].data(), Ly.RB).c_str());
    printf("  padded block 2 (rate bytes) = %s\n", hexbytes((const uint8_t*)h.h.blocks[1].data(), Ly.RB).c_str());
    printf("  reference hash recomputed from the raw bits:\n");
    printf("    final P2 input capacity   = %s, rate == A*: %s\n", HX(Ly.cap_get(h.h.p2in[29]), c),
           h.p2_ok ? "yes" : "NO");
    printf("    final S (rate, 30 words)  = %s\n", hexwords(h.h.S, 30).c_str());
    printf("    final L = %s, final R = %s\n", HX(h.h.L, c), HX(h.h.Rr, c));
    printf("    digest (first %d bits of S) = %s\n", c, hexbytes(h.h.digest.data(), h.h.digest.size()).c_str());
    if (outprefix) {
      std::string fn = std::string(outprefix) + (t ? "_Mprime.hex" : "_M.hex");
      FILE* f = fopen(fn.c_str(), "w");
      if (f) { fprintf(f, "%" PRIu64 " %s\n", h.bits, hexbytes(h.msg.data(), h.msg.size()).c_str()); fclose(f); }
    }
  }
  bool ok = true;
  auto check = [&](bool cond, const char* what) {
    printf("  [%s] %s\n", cond ? "PASS" : "FAIL", what);
    ok = ok && cond;
  };
  printf("\nAssertions:\n");
  check(H[0].msg != H[1].msg, "M != M'");
  check(A.key_ok && B.key_ok, "both joins satisfy l ^ tau_Fin = y");
  check(A.pad_ok && B.pad_ok, "both second blocks end in the padding bits 11");
  check(A.D == B.D, "D = D'");
  check(H[0].blocks_ok && H[1].blocks_ok, "reference padding gives back (m0, m1) for both messages");
  check(H[0].p1_ok && H[1].p1_ok, "final P1 inputs are u || l, so the final P1 outputs are A* || z");
  check(H[0].p2_ok && H[1].p2_ok, "final P2 inputs equal A* || D");
  check(std::memcmp(H[0].h.p2in, H[1].h.p2in, 240) == 0, "the two final P2 inputs are identical (1920 bits)");
  check(std::memcmp(H[0].h.S, H[1].h.S, 240) == 0, "final rate words S are identical (all 1920-c bits)");
  check(H[0].h.L == H[1].h.L, "final L are identical");
  check((H[0].h.Rr ^ H[1].h.Rr) == ((A.z ^ B.z) & Ly.cmask), "R ^ R' = z ^ z' (R = z ^ C differs, as expected)");
  check(H[0].h.digest == H[1].h.digest, "digests are identical");
  printf("\n%s\n", ok ? "ASSERTION PASSED: TaiChi-c(M) == TaiChi-c(M') with M != M', real P1920."
                      : "ASSERTION FAILED");
  return ok;
}

static bool verify_survivors(const Instance& I, const std::vector<Surv>& S, size_t limit) {
  const size_t n = std::min(limit, S.size());
  size_t good = 0;
  for (size_t k = 0; k < n; k++) {
    const JoinInfo J = explain_join(I, S[k].i, S[k].j);
    const HashCheck H = hash_check(I, J);
    if (J.key_ok && J.pad_ok && J.D == S[k].D && H.blocks_ok && H.p1_ok && H.p2_ok) good++;
  }
  printf("Survivor check: %zu of %zu checked survivors reproduce u || l and A* || D in the reference hash"
         " (%zu survivors in total)\n", good, n, S.size());
  return good == n;
}

// =========================================================================================
// Theory: joins q_f q_i / 2^c, survivors / 4, lambda = N^2 / 2^(c+1), N = q_f q_i / 2^(c+2).
// =========================================================================================
static const double EULER_GAMMA = 0.57721566490153286;
static double pred_median_log2q(int c) { return 0.75 * c + (5.0 + std::log2(std::log(2.0))) / 4.0; }
static double pred_mean_log2q(int c) { return 0.75 * c + (5.0 - EULER_GAMMA / std::log(2.0)) / 4.0; }
static double pred_lambda(double qf, double qi, int c) {
  const double N = qf * qi / std::ldexp(1.0, c + 2);
  return N * N / std::ldexp(1.0, c + 1);
}

// =========================================================================================
// Engine 1: interleaved (step k = forward root k + inverse query k), symmetric hash join,
// stop at the first D-collision. c <= 32.
// =========================================================================================
struct GrowTable {
  std::vector<u64> t;
  int lg = 0, c = 0;
  size_t cnt = 0, maxcnt = 0, msk = 0;
  void init(int c_, int lg0) {
    c = c_; lg = lg0; t.assign(size_t(1) << lg, 0); cnt = 0;
    msk = (size_t(1) << lg) - 1; maxcnt = (size_t(3) << lg) / 4;
  }
  inline size_t home(u64 key) const { return (size_t)((key << (64 - c)) >> (64 - lg)); }
  inline const u64* ptr(u64 key) const { return &t[home(key)]; }
  inline void put(u64 e) { size_t p = home(e >> 32); while (t[p]) p = (p + 1) & msk; t[p] = e; }
  inline void insert(u64 key, u64 idx) { put((key << 32) | (idx + 1)); cnt++; }
  template <class F> inline void lookup(u64 key, F&& f) const {
    size_t p = home(key);
    for (;;) {
      const u64 e = t[p];
      if (!e) return;
      if ((e >> 32) == key) f((e & 0xffffffffULL) - 1);
      p = (p + 1) & msk;
    }
  }
  void reserve_for(size_t extra) {
    while (cnt + extra > maxcnt) {
      std::vector<u64> old;
      old.swap(t);
      lg++; t.assign(size_t(1) << lg, 0);
      msk = (size_t(1) << lg) - 1; maxcnt = (size_t(3) << lg) / 4;
      for (u64 e : old) if (e) put(e);
    }
  }
};

struct IncResult {
  bool found = false;
  u64 q = 0, fP1 = 0, fP2 = 0, iP1 = 0, reevals = 0, joins = 0, surv = 0;
  u64 i1 = 0, j1 = 0, i2 = 0, j2 = 0;
  double secs = 0;
  std::vector<Surv> survivors;
};

static u64 qmax_for(int c) {
  double lim = std::ldexp(1.0, (int)std::ceil(0.75 * c + 4));
  lim = std::min(lim, std::ldexp(1.0, c));
  lim = std::min(lim, 4294967000.0);
  return (u64)lim;
}

static IncResult run_incremental(const Instance& I, u64 qmax, bool keep_surv) {
  const Layout& Ly = I.Ly;
  IncResult R;
  GrowTable TF, TI;
  TF.init(Ly.c, 12);
  TI.init(Ly.c, 12);
  std::unordered_map<u64, std::pair<u64, u64>> TD;
  alignas(64) u64 fx[BS], fy[BS], fsw[BS], ikey[BS], iuw[BS];
  bool found = false;
  const double t0 = now_sec();
  auto join = [&](u64 i, u64 j, u64 x, u64 y, u64 sw, u64 uw, u64 z) {
    R.joins++;
    if (((sw ^ uw) & Ly.padmask) != Ly.padmask) return;
    R.surv++;
    const u64 D = (x ^ y ^ z) & Ly.cmask;
    if (keep_surv) R.survivors.push_back({D, i, j});
    auto it = TD.find(D);
    if (it != TD.end()) {
      R.i1 = it->second.first; R.j1 = it->second.second; R.i2 = i; R.j2 = j;
      found = true;
      return;
    }
    TD.emplace(D, std::make_pair(i, j));
  };
  u64 k0 = 0;
  while (!found && k0 < qmax) {
    const int m = (int)std::min<u64>((u64)BS, qmax - k0);
    roots_range(I, k0, m, fy, fx, fsw);
    queries_range(I, k0, m, ikey, iuw);
    TF.reserve_for(m);
    TI.reserve_for(m);
    for (int t = 0; t < m; t++) {
      __builtin_prefetch(TI.ptr(fy[t]), 0, 1);
      __builtin_prefetch(TF.ptr(fy[t]), 1, 1);
      __builtin_prefetch(TF.ptr(ikey[t]), 0, 1);
      __builtin_prefetch(TI.ptr(ikey[t]), 1, 1);
    }
    for (int t = 0; t < m && !found; t++) {
      const u64 k = k0 + t;
      R.fP1++; R.fP2++;
      TI.lookup(fy[t], [&](u64 j) {
        if (found) return;
        u64 uw;
        if (j >= k0) uw = iuw[j - k0];
        else { u64 u[30], l; I.inv_query(j, u, l); uw = u[Ly.pad_word]; R.reevals++; }
        join(k, j, fx[t], fy[t], fsw[t], uw, I.z_of(j));
      });
      if (found) { R.q = k + 1; break; }
      TF.insert(fy[t], k);
      R.iP1++;
      TI.lookup(0, [](u64) {});  // no-op keeps symmetry readable
      TF.lookup(ikey[t], [&](u64 i) {
        if (found) return;
        u64 x, y, sw;
        if (i >= k0) { x = fx[i - k0]; y = fy[i - k0]; sw = fsw[i - k0]; }
        else { u64 s[30]; I.fwd_root(i, s, x, y); sw = s[Ly.pad_word]; R.reevals += 2; }
        join(i, k, x, y, sw, iuw[t], I.z_of(k));
      });
      if (found) { R.q = k + 1; break; }
      TI.insert(ikey[t], k);
    }
    k0 += m;
  }
  R.found = found;
  if (!found) R.q = k0;
  R.secs = now_sec() - t0;
  return R;
}

// =========================================================================================
// Engine 2: q_f forward roots in a compact table (8-byte slots: 32-bit fingerprint of y,
// 32-bit index), q_i inverse queries streamed against it. The y-space can be split into P
// passes when q_f does not fit in memory. Survivors go to a shared D-map, so the run can stop
// at the first D-collision. Matches are confirmed by recomputing the forward root.
// =========================================================================================
struct BigTable {
  u64* t = nullptr;
  size_t n = 0, bytes = 0;
  BigTable() = default;
  BigTable(const BigTable&) = delete;
  BigTable& operator=(const BigTable&) = delete;
  ~BigTable() { if (t) free(t); }
  void alloc(size_t nslots) {
    n = nslots;
    const size_t HP = size_t(1) << 21;
    bytes = (n * sizeof(u64) + HP - 1) / HP * HP;
    t = (u64*)aligned_alloc(HP, bytes);
    if (!t) { void* q = nullptr; if (posix_memalign(&q, 4096, bytes) == 0) t = (u64*)q; }
    if (!t) t = (u64*)malloc(bytes);
    if (!t) { fprintf(stderr, "cannot allocate %.2f GB for the join table; lower --mem-gb\n", bytes / 1e9); exit(1); }
#ifdef MADV_HUGEPAGE
    madvise(t, bytes, MADV_HUGEPAGE);
#endif
  }
  void clear(int threads) {
    std::vector<std::thread> th;
    const size_t chunk = ((bytes / threads) + 4095) & ~size_t(4095);
    for (int k = 0; k < threads; k++)
      th.emplace_back([this, k, chunk] {
        const size_t lo = std::min(bytes, (size_t)k * chunk), hi = std::min(bytes, lo + chunk);
        if (lo < hi) std::memset((char*)t + lo, 0, hi - lo);
      });
    for (auto& x : th) x.join();
  }
  inline size_t home(u64 frac) const { return (size_t)(((u128)frac * n) >> 64); }
  inline void insert(size_t p, u64 e) {
    for (;;) {
      const u64 cur = __atomic_load_n(&t[p], __ATOMIC_RELAXED);
      if (cur == 0) {
        u64 expected = 0;
        if (__atomic_compare_exchange_n(&t[p], &expected, e, false, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return;
        continue;
      }
      if (++p == n) p = 0;
    }
  }
  template <class F> inline void lookup(size_t p, u32 fp, F&& f) const {
    for (;;) {
      const u64 cur = t[p];
      if (cur == 0) return;
      if ((u32)(cur >> 32) == fp) f((u32)cur - 1);
      if (++p == n) p = 0;
    }
  }
};

struct BatchResult {
  u64 Qf = 0, Qi = 0;
  int P = 1;
  size_t nslots = 0;
  u64 fP1 = 0, fP2 = 0;          // distinct forward queries
  u64 iP1 = 0;                   // inverse queries actually made (less than Qi on early stop)
  u64 reeval_pass = 0, reeval_join = 0;
  u64 cands = 0, falsepos = 0, joins = 0, surv = 0;
  bool found = false;
  u64 i1 = 0, j1 = 0, i2 = 0, j2 = 0, Dcol = 0;
  std::vector<Surv> survivors;
  double secs = 0, secs_fwd = 0;
};

struct Cand { u64 i, j, key, uw; };

static BatchResult run_batch(const Instance& I, u64 Qf, u64 Qi, int P, double alpha, int threads, bool stop_first,
                             bool verbose) {
  const Layout& Ly = I.Ly;
  const int c = Ly.c;
  BatchResult B;
  B.Qf = Qf; B.Qi = Qi; B.P = P;
  const double per_pass = (double)Qf / P;
  B.nslots = (size_t)std::ceil(per_pass / alpha + 6.0 * std::sqrt(per_pass) + 4096.0);
  BigTable Tb;
  Tb.alloc(B.nslots);
  // Striped D-map so survivor handling does not serialise across many threads. Equal D always
  // hashes to the same stripe, so collision detection stays correct. Per-thread survivor
  // buffers are merged at the end (no lock on the hot path except the rare survivor insert).
  static const int NSTRIPE = 256;
  std::mutex mu[NSTRIPE];
  std::unordered_map<u64, std::pair<u64, u64>> TD[NSTRIPE];
  std::mutex result_mu;
  std::vector<std::vector<Surv>> local_sv(threads);
  std::atomic<bool> stop{false};
  std::atomic<u64> inv_done{0}, joins{0}, surv{0}, cands{0}, falsepos{0};
  const double t0 = now_sec();
  if (verbose)
    printf("Batch engine: q_f = %" PRIu64 " (2^%.2f), q_i = %" PRIu64 " (2^%.2f), %d pass(es), table %zu slots"
           " (%.2f GB), %d thread(s), %s\n",
           Qf, std::log2((double)Qf), Qi, std::log2((double)Qi), P, B.nslots, B.nslots * 8.0 / 1e9, threads,
           stop_first ? "stop at first D-collision" : "run to completion");
  for (int p = 0; p < P && !stop.load(); p++) {
    const double tp = now_sec();
    Tb.clear(threads);
    auto fwd_worker = [&](int k) {
      const u64 lo = Qf * (u64)k / threads, hi = Qf * (u64)(k + 1) / threads;
      alignas(64) u64 y[BS];
      size_t hs[BS];
      u64 es[BS];
      for (u64 i0 = lo; i0 < hi; i0 += BS) {
        const int m = (int)std::min<u64>((u64)BS, hi - i0);
        roots_range(I, i0, m, y, nullptr, nullptr);
        int cnt = 0;
        for (int t = 0; t < m; t++) {
          const u128 tt = (u128)(y[t] << (64 - c)) * (u64)P;
          if ((int)(u64)(tt >> 64) != p) continue;
          const size_t h = Tb.home((u64)tt);
          __builtin_prefetch(&Tb.t[h], 1, 0);
          hs[cnt] = h;
          es[cnt] = ((u64)(u32)y[t] << 32) | (i0 + t + 1);
          cnt++;
        }
        for (int q = 0; q < cnt; q++) Tb.insert(hs[q], es[q]);
      }
    };
    {
      std::vector<std::thread> th;
      for (int k = 0; k < threads; k++) th.emplace_back(fwd_worker, k);
      for (auto& x : th) x.join();
    }
    const double tf = now_sec();
    B.secs_fwd += tf - tp;
    if (verbose) {
      printf("  pass %d/%d: forward table built in %.1f s\n", p + 1, P, tf - tp);
      fflush(stdout);
    }
    auto inv_worker = [&](int k) {
      const u64 lo = Qi * (u64)k / threads, hi = Qi * (u64)(k + 1) / threads;
      alignas(64) u64 key[BS], uw[BS];
      alignas(64) u64 cidx[BS], cy[BS], cx[BS], csw[BS];
      size_t hs[BS];
      int ts[BS];
      std::vector<Cand> cand;
      cand.reserve(4 * BS);
      std::vector<Surv>& local = local_sv[k];
      double last_report = now_sec();
      auto flush = [&]() {
        for (size_t done = 0; done < cand.size();) {
          const int m = (int)std::min<size_t>((size_t)BS, cand.size() - done);
          for (int t = 0; t < m; t++) cidx[t] = cand[done + t].i;
          roots_idx(I, cidx, m, cy, cx, csw);
          for (int t = 0; t < m; t++) {
            const Cand& C = cand[done + t];
            cands++;
            if (cy[t] != C.key) { falsepos++; continue; }
            joins++;
            if (((csw[t] ^ C.uw) & Ly.padmask) != Ly.padmask) continue;
            surv++;
            const u64 D = (cx[t] ^ cy[t] ^ I.z_of(C.j)) & Ly.cmask;
            local.push_back({D, C.i, C.j});
            const int st = (int)(D & (NSTRIPE - 1));
            std::lock_guard<std::mutex> g(mu[st]);
            auto it = TD[st].find(D);
            if (it != TD[st].end() && !(it->second.first == C.i && it->second.second == C.j)) {
              std::lock_guard<std::mutex> rg(result_mu);
              if (!B.found) {
                B.found = true;
                B.i1 = it->second.first; B.j1 = it->second.second; B.i2 = C.i; B.j2 = C.j; B.Dcol = D;
              }
              if (stop_first) stop.store(true);
            } else {
              TD[st].emplace(D, std::make_pair(C.i, C.j));
            }
          }
          done += m;
        }
        cand.clear();
      };
      for (u64 j0 = lo; j0 < hi && !stop.load(std::memory_order_relaxed); j0 += BS) {
        const int m = (int)std::min<u64>((u64)BS, hi - j0);
        queries_range(I, j0, m, key, uw);
        int cnt = 0;
        for (int t = 0; t < m; t++) {
          const u128 tt = (u128)(key[t] << (64 - c)) * (u64)P;
          if ((int)(u64)(tt >> 64) != p) continue;
          const size_t h = Tb.home((u64)tt);
          __builtin_prefetch(&Tb.t[h], 0, 0);
          hs[cnt] = h;
          ts[cnt] = t;
          cnt++;
        }
        for (int q = 0; q < cnt; q++) {
          const int t = ts[q];
          Tb.lookup(hs[q], (u32)key[t], [&](u32 i) { cand.push_back({(u64)i, j0 + t, key[t], uw[t]}); });
        }
        if (cand.size() >= (size_t)BS) flush();
        if (p == 0) inv_done.fetch_add((u64)m, std::memory_order_relaxed);
        if (verbose && k == 0 && now_sec() - last_report > 60) {
          last_report = now_sec();
          const double frac = (double)inv_done.load() / (double)Qi;
          const double el = last_report - tf;
          printf("    inverse queries %5.1f%% (%" PRIu64 "), joins %" PRIu64 ", survivors %" PRIu64
                 ", elapsed %.0f s, full pass ETA %.0f s\n",
                 100 * frac, inv_done.load(), joins.load(), surv.load(), el, frac > 0 ? el / frac - el : 0.0);
          fflush(stdout);
        }
      }
      flush();
    };
    {
      std::vector<std::thread> th;
      for (int k = 0; k < threads; k++) th.emplace_back(inv_worker, k);
      for (auto& x : th) x.join();
    }
    if (verbose) {
      printf("  pass %d/%d: inverse phase %.1f s, joins %" PRIu64 ", survivors %" PRIu64 "%s\n", p + 1, P,
             now_sec() - tf, joins.load(), surv.load(), B.found ? ", D-collision found" : "");
      fflush(stdout);
    }
  }
  for (auto& v : local_sv) B.survivors.insert(B.survivors.end(), v.begin(), v.end());
  B.cands = cands; B.falsepos = falsepos; B.joins = joins; B.surv = surv;
  B.fP1 = Qf; B.fP2 = Qf;
  B.iP1 = (P == 1) ? inv_done.load() : Qi;
  B.reeval_pass = (u64)(P - 1) * (2 * Qf + Qi);
  B.reeval_join = 2 * B.cands;
  B.secs = now_sec() - t0;
  return B;
}

struct ColPair { u64 i1, j1, i2, j2, D, qint; };
static std::vector<ColPair> find_collisions(std::vector<Surv> S) {
  std::sort(S.begin(), S.end(), [](const Surv& a, const Surv& b) {
    if (a.D != b.D) return a.D < b.D;
    if (a.i != b.i) return a.i < b.i;
    return a.j < b.j;
  });
  std::vector<ColPair> out;
  for (size_t a = 0; a < S.size();) {
    size_t b = a + 1;
    while (b < S.size() && S[b].D == S[a].D) b++;
    for (size_t x = a; x < b; x++)
      for (size_t y = x + 1; y < b; y++) {
        ColPair cp{S[x].i, S[x].j, S[y].i, S[y].j, S[x].D, 0};
        cp.qint = std::max(std::max(cp.i1, cp.j1), std::max(cp.i2, cp.j2)) + 1;
        out.push_back(cp);
      }
    a = b;
  }
  std::sort(out.begin(), out.end(), [](const ColPair& a, const ColPair& b) { return a.qint < b.qint; });
  return out;
}

// =========================================================================================
// Command line
// =========================================================================================
struct Args {
  std::string mode;
  int c = 40, threads = 0, trials = 0, passes = 0;   // threads 0 => all hardware cores
  size_t verify = 2000;
  u64 seed = 0x7a1c41a77ac4c0deULL;
  double mem_gb = 5.0, alpha = 0.75, lambda = 3.0, logq = -1;
  u64 q = 0, qf = 0, qi = 0;
  bool stop_first = true;
  std::string cs = "8,16,24,32";
  std::string csv, out;
};
static void usage() {
  printf("usage: taichi_real <selftest|bench|single|scale|batch> [options]\n"
         "  --c C (multiple of 8, 8..64)  --seed S  --threads T\n"
         "  scale : --cs 8,16,24,32 --trials N --csv FILE\n"
         "  batch : --lambda L (target expected D-collisions, default 3) --mem-gb G (table, default 5)\n"
         "          --q Q | --logq L (balanced)  or  --qf Q --qi Q   --passes P  --no-stop  --out PREFIX\n");
}
static Args parse_args(int argc, char** argv) {
  Args A;
  if (argc < 2) { usage(); exit(1); }
  A.mode = argv[1];
  for (int k = 2; k < argc; k++) {
    std::string o = argv[k];
    auto val = [&]() -> std::string {
      if (k + 1 >= argc) { fprintf(stderr, "missing value for %s\n", o.c_str()); exit(1); }
      return argv[++k];
    };
    if (o == "--c") A.c = std::stoi(val());
    else if (o == "--threads") A.threads = std::stoi(val());
    else if (o == "--trials") A.trials = std::stoi(val());
    else if (o == "--passes") A.passes = std::stoi(val());
    else if (o == "--verify") A.verify = (size_t)std::stoull(val());
    else if (o == "--seed") A.seed = std::stoull(val(), nullptr, 0);
    else if (o == "--mem-gb") A.mem_gb = std::stod(val());
    else if (o == "--alpha") A.alpha = std::stod(val());
    else if (o == "--lambda") A.lambda = std::stod(val());
    else if (o == "--logq") A.logq = std::stod(val());
    else if (o == "--q") A.q = std::stoull(val());
    else if (o == "--qf") A.qf = std::stoull(val());
    else if (o == "--qi") A.qi = std::stoull(val());
    else if (o == "--no-stop") A.stop_first = false;
    else if (o == "--cs") A.cs = val();
    else if (o == "--csv") A.csv = val();
    else if (o == "--out") A.out = val();
    else { fprintf(stderr, "unknown option %s\n", o.c_str()); usage(); exit(1); }
  }
  if (A.threads <= 0) {
    unsigned hc = std::thread::hardware_concurrency();
    A.threads = hc ? (int)hc : 4;
  }
  return A;
}
static void check_c(int c) {
  if (c < 8 || c > 64 || c % 8) { fprintf(stderr, "c must be a multiple of 8 in [8, 64]\n"); exit(1); }
}
// Measured forward / inverse P1920 throughput (calls/s) using `threads` threads, for ETAs.
static void calibrate(const Instance& I, int threads, double& fwd_rate, double& inv_rate) {
  auto probe = [&](bool inv) {
    std::atomic<u64> total{0};
    std::vector<std::thread> th;
    const double secs = 0.4;
    for (int k = 0; k < threads; k++)
      th.emplace_back([&, k] {
        alignas(64) u64 st[30][BS];
        SplitMix64 g(0x1234 + k);
        for (int w = 0; w < 30; w++) for (int l = 0; l < BS; l++) st[w][l] = g.next();
        u64 n = 0;
        const double t0 = now_sec();
        while (now_sec() - t0 < secs) {
          for (int r = 0; r < 8; r++) { if (inv) P1920_inv_batch(st, 1); else P1920_fwd_batch(st, 1); }
          n += 8 * BS;
        }
        total += n;
      });
    double t0 = now_sec();
    for (auto& x : th) x.join();
    return (double)total.load() / (now_sec() - t0);
  };
  fwd_rate = probe(false);
  inv_rate = probe(true);
}
static void print_header(const Instance& I, const char* mode) {
  printf("TaiChi IFS structural collision attack on the real P1920 | mode: %s\n", mode);
  printf("  state 1920 bits, rate r = 1920 - c = %" PRIu64 " bits (%" PRIu64 " bytes), capacity c = %d bits"
         " (high %d bits of word 29)\n", I.Ly.R, I.Ly.RB, I.Ly.c, I.Ly.c);
  printf("  P1 = P1920(., 1), P2 = P1920(., 2) from the NGCC reference; tau_Abs = 1, tau_Fin = 2; 10*1 padding;"
         " digest = first %d bits of S\n", I.Ly.c);
  printf("  seed = 0x%016" PRIx64 " (selects A*, the m0 family and the z family; the permutation is fixed)\n", I.seed);
}

// =========================================================================================
// selftest
// =========================================================================================
static int mode_selftest(const Args& A) {
  int fails = 0, printed = 0;
  auto CHECK = [&](bool cond, const char* what) {
    if (!cond) { fails++; if (printed++ < 20) printf("  FAIL: %s\n", what); }
  };
  // [1] P1920 forward against reference vectors (state = splitmix64(seed=1), banks 1 and 2)
  {
    u64 st[30], z = 1;
    for (int i = 0; i < 30; i++) { z += 0x9e3779b97f4a7c15ULL; st[i] = mix64(z); }
    u64 s1[30], s2[30];
    std::memcpy(s1, st, sizeof st);
    std::memcpy(s2, st, sizeof st);
    P1920_fwd(s1, 1);
    P1920_fwd(s2, 2);
    CHECK(s1[0] == 0xd85d46fd850dfe30ULL && s1[1] == 0xc547d13764dbbf44ULL && s1[29] == 0x2d6fb8c4e5c38104ULL,
          "P1 vector");
    CHECK(s2[0] == 0xe73802569f47c69fULL && s2[1] == 0x1030ef412b3fbed1ULL && s2[29] == 0x0f0d1f9840f53eb9ULL,
          "P2 vector");
  }
  printf("[1] P1920 forward matches the reference vectors (both banks)\n");
  // [2] real TaiChi-512 digests
  {
    const char* expect[4] = {
        "3cbda17667d0ed6c468e7a4e6f732f97bd3ed7a0acde4f9831e107af120cf74f82112b3362f8eb75ac0aa0d7748ffbe080b57164a5a64f7c8bbdd7c4e3f2de34",
        "83bd7a106e5383d4900dae4868ea5731bb2d6cae0fdc0df25e12e84aa0652565c193872c44845a7a92b156177b5b500c616a43ec9acd25c43f97e757ee03e010",
        "c505510eb0ef8f1ed98ecac0cccfb36913633324581ff20ff405e1801c02f2ba6c834a19ebb64739e538382d77d74b260671e974f82cdd9cbc019b3cf10a1eb0",
        "20cff92aabe3f6b40b5f5f14e2ea1af612cd2f1ea5f93e493cb0b426ea1992096e1c8b8300f5094036b5db9071820fc35a248cc34316c22d2f88455748521ab3"};
    const char* msgs[3] = {"", "abc", "The quick brown fox jumps over the lazy dog"};
    for (int t = 0; t < 4; t++) {
      std::vector<uint8_t> m;
      if (t < 3) m.assign(msgs[t], msgs[t] + strlen(msgs[t]));
      else for (int i = 0; i < 200; i++) m.push_back((uint8_t)(i * 37 + 11));
      m.push_back(0);
      uint8_t d[64];
      real_taichi512(m.data(), (u64)(m.size() - 1) * 8, d);
      CHECK(hexbytes(d, 64) == expect[t], "real TaiChi-512 digest");
    }
  }
  printf("[2] real TaiChi-512 (rate 1408, capacity 512) reproduces the 4 reference digests\n");
  // [3] inverse
  {
    SplitMix64 g(99);
    int bad = 0;
    for (int k = 0; k < 20000; k++) {
      u64 st[30], cp[30];
      for (int i = 0; i < 30; i++) cp[i] = st[i] = g.next();
      const int b = 1 + (k & 1);
      P1920_fwd(st, b);
      P1920_inv(st, b);
      if (std::memcmp(st, cp, sizeof st)) bad++;
    }
    alignas(64) u64 Bt[30][BS], Cp[30][BS];
    for (int w = 0; w < 30; w++)
      for (int l = 0; l < BS; l++) Cp[w][l] = Bt[w][l] = g.next();
    P1920_inv_batch(Bt, 1);
    for (int l = 0; l < BS; l++) {
      u64 s[30];
      for (int w = 0; w < 30; w++) s[w] = Bt[w][l];
      P1920_fwd(s, 1);
      for (int w = 0; w < 30; w++) if (s[w] != Cp[w][l]) bad++;
    }
    CHECK(bad == 0, "inverse round trip");
  }
  printf("[3] P1920 inverse: scalar and batched round trips\n");
  for (int c : {8, 16, 24, 32, 40, 48, 64}) {
    Instance I;
    I.init(c, 0x5e1f7e57ULL + c);
    const Layout& Ly = I.Ly;
    // [4] first-block state via the reduced reference step = (s, y, x ^ y)
    for (u64 i = 0; i < 50; i++) {
      u64 s[30], x, y;
      I.fwd_root(i, s, x, y);
      u64 m0[30];
      I.m0_of(i, m0);
      // one absorb step through a 1-block message prefix: compare by running two blocks and
      // checking the block-1 chaining value through explain_join's own recomputation
      u64 st[30];
      std::memcpy(st, m0, sizeof st);
      st[29] = Ly.cap_set(st[29], 0 ^ TAU_ABS);
      P1920_fwd(st, 1);
      const u64 B = Ly.cap_get(st[29]);
      st[29] = Ly.cap_set(st[29], 0 ^ B);
      P1920_fwd(st, 2);
      const u64 C = Ly.cap_get(st[29]);
      st[29] &= Ly.lowmask;
      CHECK(B == x && C == y && std::memcmp(st, s, sizeof st) == 0, "first-block state");
    }
    // [5] padding round trip + attack equations through the reduced reference hash
    for (u64 trial = 0; trial < 20; trial++) {
      JoinInfo J = explain_join(I, trial, trial * 7 + 3);
      J.m1[Ly.pad_word] |= Ly.padmask;   // force a valid final block for the round trip
      u64 bits;
      std::vector<uint8_t> raw = raw_message(Ly, J.m0, J.m1, bits);
      RTrace h = reduced_hash(Ly, raw.data(), bits);
      CHECK(h.blocks.size() == 2 && std::memcmp(h.blocks[0].data(), J.m0, 240) == 0 &&
                std::memcmp(h.blocks[1].data(), J.m1, 240) == 0,
            "padding round trip");
    }
  }
  printf("[4] reduced mode, c in {8,16,24,32,40,48,64}: first-block state is (s, y, x^y)\n");
  printf("[5] reduced mode: 2R-2 bit messages pad back to (m0, m1) with the '11' in word %s\n", "pad_word");
  // [6] engines: partitioned batch engine == brute force, interleaved first-collision q agrees
  {
    int with_col = 0;
    for (int trial = 0; trial < 6; trial++) {
      Instance I;
      I.init(16, 777 + trial);
      const Layout& Ly = I.Ly;
      const u64 Q = 12000;
      std::vector<u64> X(Q), Y(Q), SW(Q), KEY(Q), UW(Q);
      for (u64 i = 0; i < Q; i++) { u64 s[30]; I.fwd_root(i, s, X[i], Y[i]); SW[i] = s[Ly.pad_word]; }
      for (u64 j = 0; j < Q; j++) { u64 u[30], l; I.inv_query(j, u, l); KEY[j] = l ^ TAU_FIN; UW[j] = u[Ly.pad_word]; }
      std::unordered_multimap<u64, u64> byY;
      for (u64 i = 0; i < Q; i++) byY.emplace(Y[i], i);
      std::vector<Surv> bf;
      u64 bj = 0;
      for (u64 j = 0; j < Q; j++) {
        auto rg = byY.equal_range(KEY[j]);
        for (auto it = rg.first; it != rg.second; ++it) {
          const u64 i = it->second;
          bj++;
          if (((SW[i] ^ UW[j]) & Ly.padmask) == Ly.padmask) bf.push_back({(X[i] ^ Y[i] ^ I.z_of(j)) & Ly.cmask, i, j});
        }
      }
      BatchResult B = run_batch(I, Q, Q, 3, 0.75, 2, false, false);
      CHECK(B.joins == bj, "batch joins == brute force");
      auto key = [](const Surv& a, const Surv& b) {
        if (a.D != b.D) return a.D < b.D;
        if (a.i != b.i) return a.i < b.i;
        return a.j < b.j;
      };
      std::vector<Surv> bs = B.survivors;
      std::sort(bs.begin(), bs.end(), key);
      std::sort(bf.begin(), bf.end(), key);
      bool same = bs.size() == bf.size();
      for (size_t t = 0; same && t < bs.size(); t++) same = bs[t].D == bf[t].D && bs[t].i == bf[t].i && bs[t].j == bf[t].j;
      CHECK(same, "batch survivors == brute force");
      auto cols = find_collisions(B.survivors);
      IncResult R = run_incremental(I, Q, false);
      if (R.found) {
        with_col++;
        CHECK(!cols.empty() && cols[0].qint == R.q, "first-collision q: interleaved == partitioned");
      } else {
        CHECK(cols.empty(), "no collision in either engine");
      }
    }
    printf("[6] c=16: 3-pass engine == brute force on 6 keys; first-collision q agrees (%d keys with a collision)\n",
           with_col);
  }
  // [7] a full collision at c=16 through the reduced reference hash; all survivors reproduce
  {
    Instance I;
    I.init(16, 31337);
    IncResult R = run_incremental(I, qmax_for(16), true);
    CHECK(R.found, "collision found at c=16");
    size_t good = 0;
    for (auto& sv : R.survivors) {
      const JoinInfo J = explain_join(I, sv.i, sv.j);
      const HashCheck H = hash_check(I, J);
      if (J.key_ok && J.pad_ok && J.D == sv.D && H.blocks_ok && H.p1_ok && H.p2_ok) good++;
    }
    CHECK(good == R.survivors.size(), "survivors reproduce");
    const JoinInfo J1 = explain_join(I, R.i1, R.j1), J2 = explain_join(I, R.i2, R.j2);
    const HashCheck H1 = hash_check(I, J1), H2 = hash_check(I, J2);
    CHECK(H1.msg != H2.msg && H1.h.digest == H2.h.digest && std::memcmp(H1.h.S, H2.h.S, 240) == 0,
          "c=16 collision through the reference hash");
    printf("[7] c=16: %zu/%zu survivors reproduce A*||D in the reference hash; collision at q = %" PRIu64 " verified\n",
           good, R.survivors.size(), R.q);
  }
  printf("selftest: %s (%d failures)\n", fails ? "FAILED" : "all checks passed", fails);
  return fails ? 1 : 0;
}

// =========================================================================================
// bench
// =========================================================================================
static int mode_bench(const Args& A) {
  check_c(A.c);
  Instance I;
  I.init(A.c, A.seed);
  auto run = [&](int threads, double secs, u64& roots, u64& queries) {
    std::atomic<bool> go{true};
    std::atomic<u64> r{0}, q{0};
    std::vector<std::thread> th;
    for (int k = 0; k < threads; k++)
      th.emplace_back([&, k] {
        alignas(64) u64 y[BS], key[BS], uw[BS];
        u64 base = (u64)k << 40;
        while (go.load(std::memory_order_relaxed)) {
          roots_range(I, base, BS, y, nullptr, nullptr);
          queries_range(I, base, BS, key, uw);
          base += BS;
          r += BS;
          q += BS;
        }
      });
    std::this_thread::sleep_for(std::chrono::duration<double>(secs));
    go = false;
    for (auto& x : th) x.join();
    roots = r;
    queries = q;
  };
  for (int t : {1, A.threads}) {
    u64 r, q;
    const double secs = 4.0;
    run(t, secs, r, q);
    printf("threads=%d: %.3e (root + query) steps/s, i.e. %.0f ns per step (2 forward + 1 inverse P1920)\n", t,
           r / secs, secs * 1e9 / r);
  }
  return 0;
}

// =========================================================================================
// single: interleaved run with the full report (c <= 32)
// =========================================================================================
static int mode_single(const Args& A) {
  check_c(A.c);
  if (A.c > 32) { fprintf(stderr, "single mode supports c <= 32; use batch\n"); return 1; }
  Instance I;
  I.init(A.c, A.seed);
  print_header(I, "single (interleaved, stop at the first D-collision)");
  IncResult R = run_incremental(I, qmax_for(A.c), true);
  const int c = A.c;
  const double q = (double)R.q;
  printf("\nResults:\n");
  printf("  r = %" PRIu64 ", c = %d, q = %" PRIu64 " (log2 q = %.3f, log2 q / c = %.4f)\n", I.Ly.R, c, R.q,
         std::log2(q), std::log2(q) / c);
  printf("  forward P1 queries (first block)       : %" PRIu64 "\n", R.fP1);
  printf("  first-stage P2 queries                 : %" PRIu64 "\n", R.fP2);
  printf("  inverse P1 queries                     : %" PRIu64 "\n", R.iP1);
  printf("  repeated evaluations for join recovery : %" PRIu64 "\n", R.reevals);
  printf("  valid forward/inverse joins            : %" PRIu64 " (predicted %.1f)\n", R.joins,
         (double)R.fP1 * (double)R.iP1 / std::ldexp(1.0, c));
  printf("  joins surviving padding                : %" PRIu64 " (predicted %.1f)\n", R.surv,
         (double)R.fP1 * (double)R.iP1 / std::ldexp(1.0, c + 2));
  printf("  predicted median first-collision q     : 2^%.3f; this run 2^%.3f\n", pred_median_log2q(c), std::log2(q));
  printf("  time                                   : %.2f s\n", R.secs);
  if (!R.found) { printf("no collision\n"); return 1; }
  const JoinInfo J1 = explain_join(I, R.i1, R.j1), J2 = explain_join(I, R.i2, R.j2);
  const bool ok = report_collision(I, J1, J2, A.out.empty() ? nullptr : A.out.c_str());
  printf("\n");
  const bool sok = verify_survivors(I, R.survivors, A.verify);
  return (ok && sok) ? 0 : 2;
}

// =========================================================================================
// scale
// =========================================================================================
static int default_trials(int c) {
  if (c <= 8) return 2000;
  if (c <= 16) return 1000;
  if (c <= 24) return 300;
  return 24;
}
static int mode_scale(const Args& A) {
  std::vector<int> cs;
  for (size_t p = 0; p < A.cs.size();) {
    size_t e = A.cs.find(',', p);
    if (e == std::string::npos) e = A.cs.size();
    cs.push_back(std::stoi(A.cs.substr(p, e - p)));
    p = e + 1;
  }
  printf("Scaling of the first-collision q on the real P1920 (interleaved: step k = root k + query k)\n");
  printf("  each trial: fresh A*, m0 family and z family (the permutation is the real one, fixed)\n");
  printf("  theory: median log2 q = 0.75c + %.4f, mean = 0.75c + %.4f, sd %.3f\n\n", pred_median_log2q(0),
         pred_mean_log2q(0), (M_PI / std::sqrt(6.0)) / std::log(2.0) / 4.0);
  printf("   c  trials  med log2q  mean log2q    se    pred med  pred mean  log2(q_med)/c  P[q<=pred med]  joins/pred  surv/joins   time\n");
  std::vector<double> X, Y, SE;
  FILE* csv = A.csv.empty() ? nullptr : fopen(A.csv.c_str(), "w");
  if (csv) fprintf(csv, "c,trials,median_log2q,mean_log2q,se,sd,pred_median,pred_mean,frac_le_pred_median,joins_ratio,surv_ratio,seconds\n");
  for (int c : cs) {
    check_c(c);
    if (c > 32) { printf("  c = %d skipped (interleaved engine is for c <= 32)\n", c); continue; }
    const int T = A.trials > 0 ? A.trials : default_trials(c);
    std::vector<IncResult> res(T);
    std::atomic<int> next{0};
    const double t0 = now_sec();
    auto worker = [&]() {
      for (;;) {
        const int t = next.fetch_add(1);
        if (t >= T) return;
        Instance I;
        I.init(c, mix64(A.seed ^ mix64(((u64)c << 40) ^ (u64)t)));
        res[t] = run_incremental(I, qmax_for(c), false);
      }
    };
    std::vector<std::thread> th;
    for (int k = 0; k < A.threads; k++) th.emplace_back(worker);
    for (auto& x : th) x.join();
    const double secs = now_sec() - t0;
    std::vector<double> lq;
    double jr = 0, sr = 0;
    for (auto& R : res) {
      if (!R.found) continue;
      lq.push_back(std::log2((double)R.q));
      jr += R.joins / ((double)R.fP1 * (double)R.iP1 / std::ldexp(1.0, c));
      sr += R.joins ? (double)R.surv / R.joins : 0;
    }
    if (lq.empty()) { printf("  c = %d: no collisions\n", c); continue; }
    std::sort(lq.begin(), lq.end());
    const size_t n = lq.size();
    const double med = (n % 2) ? lq[n / 2] : 0.5 * (lq[n / 2 - 1] + lq[n / 2]);
    double mean = 0, var = 0;
    for (double v : lq) mean += v;
    mean /= n;
    for (double v : lq) var += (v - mean) * (v - mean);
    const double sd = std::sqrt(var / std::max<size_t>(1, n - 1)), se = sd / std::sqrt((double)n);
    int le = 0;
    for (double v : lq) le += (v <= pred_median_log2q(c));
    printf("  %2d  %6zu  %9.3f  %10.3f  %6.3f  %8.3f  %9.3f  %13.4f  %14.3f  %10.3f  %10.3f  %6.1fs\n", c, n, med,
           mean, se, pred_median_log2q(c), pred_mean_log2q(c), med / c, (double)le / n, jr / n, sr / n, secs);
    fflush(stdout);
    if (csv)
      fprintf(csv, "%d,%zu,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.1f\n", c, n, med, mean, se, sd,
              pred_median_log2q(c), pred_mean_log2q(c), (double)le / n, jr / n, sr / n, secs);
    X.push_back(c); Y.push_back(mean); SE.push_back(se);
  }
  if (csv) fclose(csv);
  if (X.size() >= 2) {
    double Sw = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for (size_t k = 0; k < X.size(); k++) {
      const double w = 1.0 / (SE[k] * SE[k]);
      Sw += w; Sx += w * X[k]; Sy += w * Y[k]; Sxx += w * X[k] * X[k]; Sxy += w * X[k] * Y[k];
    }
    const double Dt = Sw * Sxx - Sx * Sx;
    const double s = (Sw * Sxy - Sx * Sy) / Dt, b = (Sxx * Sy - Sx * Sxy) / Dt;
    printf("\nWeighted fit of mean log2 q on c: slope = %.4f +- %.4f (theory 0.75), intercept = %.3f +- %.3f (theory %.3f)\n",
           s, std::sqrt(Sw / Dt), b, std::sqrt(Sxx / Dt), pred_mean_log2q(0));
    double chi = 0;
    for (size_t k = 0; k < X.size(); k++) chi += std::pow((Y[k] - pred_mean_log2q((int)X[k])) / SE[k], 2);
    printf("chi2 of the parameter-free prediction 0.75c + %.3f: %.2f over %zu points\n", pred_mean_log2q(0), chi, X.size());
  }
  return 0;
}

// =========================================================================================
// batch: large c
// =========================================================================================
static int mode_batch(const Args& A) {
  check_c(A.c);
  const int c = A.c;
  Instance I;
  I.init(c, A.seed);
  const double slots = A.mem_gb * 1e9 / 8.0;
  const double M = std::floor(slots * A.alpha);            // roots per pass that fit in memory
  // product needed for lambda expected D-collisions: q_f q_i = sqrt(lambda 2^(c+1)) 2^(c+2)
  const double K = std::sqrt(A.lambda * std::ldexp(1.0, c + 1)) * std::ldexp(1.0, c + 2);
  double qf, qi;
  std::string plan;
  if (A.q || A.logq > 0) {
    qf = qi = A.q ? (double)A.q : std::pow(2.0, A.logq);
    plan = "balanced (user)";
  } else if (A.qf && A.qi) {
    qf = (double)A.qf; qi = (double)A.qi;
    plan = "user q_f, q_i";
  } else {
    const double qb = std::sqrt(K);
    if (qb <= M) { qf = qi = qb; plan = "balanced, one pass"; }
    else { qf = M; qi = K / M; plan = "memory-limited: q_f = table capacity, q_i = K / q_f, one pass"; }
  }
  const u64 Qf = (u64)std::llround(qf), Qi = (u64)std::llround(qi);
  if (Qf >= 0xffffffffULL) { fprintf(stderr, "q_f must be < 2^32 - 1\n"); return 1; }
  if ((double)Qi > std::ldexp(1.0, c)) { fprintf(stderr, "q_i must be <= 2^c\n"); return 1; }
  const int P = A.passes > 0 ? A.passes : std::max(1, (int)std::ceil((double)Qf / M));
  print_header(I, "batch");
  const double qbal = std::sqrt(std::sqrt(A.lambda * std::ldexp(1.0, c + 1)) * std::ldexp(1.0, c + 2));
  printf("  plan: %s; target lambda = %.2f (P[success] = %.3f)\n", plan.c_str(), pred_lambda((double)Qf, (double)Qi, c),
         1 - std::exp(-pred_lambda((double)Qf, (double)Qi, c)));
  printf("  balanced optimum for this lambda: q_f = q_i = 2^%.3f, total %.3g calls (2^%.3f); this plan: total %.3g"
         " calls (2^%.3f)\n", std::log2(qbal), 3 * qbal, std::log2(3 * qbal), (double)P * (2.0 * Qf + Qi),
         std::log2((double)P * (2.0 * Qf + Qi)));
  double fr, ir;
  calibrate(I, A.threads, fr, ir);
  const double eta = P * (2.0 * Qf / fr + (double)Qi / ir);
  printf("  calibration on %d threads: forward %.1f Mcalls/s, inverse %.1f Mcalls/s; estimated wall-clock %.1f min"
         " (worst case, no early stop)\n", A.threads, fr / 1e6, ir / 1e6, eta / 60);
  fflush(stdout);
  BatchResult B = run_batch(I, Qf, Qi, P, A.alpha, A.threads, A.stop_first, true);
  const double qi_used = (double)B.iP1;
  printf("\nResults:\n");
  printf("  r = %" PRIu64 ", c = %d, q_f = %" PRIu64 " (2^%.3f), q_i planned = %" PRIu64 " (2^%.3f)\n", I.Ly.R, c, Qf,
         std::log2((double)Qf), Qi, std::log2((double)Qi));
  printf("  forward P1 queries (first block)       : %" PRIu64 "\n", B.fP1);
  printf("  first-stage P2 queries                 : %" PRIu64 "\n", B.fP2);
  printf("  inverse P1 queries                     : %" PRIu64 " (2^%.3f)%s\n", B.iP1, std::log2(qi_used),
         B.found && A.stop_first ? ", stopped at the first D-collision" : "");
  printf("  repeated evaluations                   : %" PRIu64 " for memory passes + %" PRIu64 " for join recovery\n",
         B.reeval_pass, B.reeval_join);
  printf("  table hits / fingerprint false positives: %" PRIu64 " / %" PRIu64 "\n", B.cands, B.falsepos);
  printf("  valid forward/inverse joins            : %" PRIu64 " (predicted %.1f)\n", B.joins,
         (double)Qf * qi_used / std::ldexp(1.0, c));
  printf("  joins surviving padding                : %" PRIu64 " (predicted %.1f)\n", B.surv,
         (double)Qf * qi_used / std::ldexp(1.0, c + 2));
  printf("  expected D-collisions at this point    : %.3f\n", pred_lambda((double)Qf, qi_used, c));
  printf("  time                                   : %.1f s (forward table %.1f s)\n", B.secs, B.secs_fwd);
  printf("  total permutation calls                : 2^%.3f\n", std::log2(2.0 * Qf * P + qi_used + 2.0 * B.cands));
  if (!B.found) {
    printf("\nNo D-collision. Rerun with a larger --lambda or another --seed.\n");
    verify_survivors(I, B.survivors, A.verify);
    return 1;
  }
  const JoinInfo J1 = explain_join(I, B.i1, B.j1), J2 = explain_join(I, B.i2, B.j2);
  const bool ok = report_collision(I, J1, J2, A.out.empty() ? nullptr : A.out.c_str());
  printf("\n");
  const bool sok = verify_survivors(I, B.survivors, A.verify);
  return (ok && sok) ? 0 : 2;
}

int main(int argc, char** argv) {
  const Args A = parse_args(argc, argv);
  if (A.mode == "selftest") return mode_selftest(A);
  if (A.mode == "bench") return mode_bench(A);
  if (A.mode == "single") return mode_single(A);
  if (A.mode == "scale") return mode_scale(A);
  if (A.mode == "batch") return mode_batch(A);
  usage();
  return 1;
}
