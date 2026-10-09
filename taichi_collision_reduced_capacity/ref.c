/* Reference TaiChi-512 (verbatim from the NGCC submission), plus a small driver.
   Used only as ground truth to validate the reimplementation in the attack. */
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define RATE_BITS    1408
#define DIGEST_BYTES 64
#define RATE_BYTES     (RATE_BITS / 8)
#define RATE_WORDS     (RATE_BITS / 64)
#define CAPACITY_WORDS ((1920 - RATE_BITS) / 64)

typedef enum { DOMAIN_AB = 1, DOMAIN_FIN = 2, DOMAIN_SQ = 3 } DOMAIN;

static inline uint64_t ROTL64(uint64_t x, unsigned int n) { return (x << n) | (x >> ((64 - n) & 63)); }

static uint64_t RC[120];
static void precompute_RC() {
    static const int p[7] = {0, 2, 3, 6, 13, 28, 59};
    uint8_t s[7] = {1, 0, 0, 0, 0, 0, 0};
    for (int q = 0; q < 120; q++) {
        uint64_t C = 0;
        for (int i = 0; i < 7; i++) if (s[i]) C ^= (1ULL << p[i]);
        RC[q] = C;
        uint8_t next_s[7];
        next_s[0] = s[6];
        next_s[1] = s[0] ^ s[6];
        next_s[2] = s[1];
        next_s[3] = s[2];
        next_s[4] = s[3];
        next_s[5] = s[4];
        next_s[6] = s[5];
        for (int i = 0; i < 7; i++) s[i] = next_s[i];
    }
}
static void P1920(uint64_t A[30], int b) {
    static const uint8_t n1[16] = {0, 1, 2, 4, 3, 5, 6, 7, 13, 14, 15, 16, 15, 18, 11, 20};
    static const uint8_t n2[16] = {1, 2, 4, 3, 5, 6, 0, 11, 8, 9, 10, 12, 17, 14, 19, 13};
    static const uint8_t m[16] = {8, 0, 0, 0, 0, 0, 23, 62, 35, 14, 48, 1, 57, 63, 58, 22};
    uint64_t tmp[30];
    uint64_t *cur = A;
    uint64_t *next = tmp;
    int t;
    precompute_RC();
    for (t = 0; t < 12; t++) {
        int j;
        for (j = 0; j < 5; j++) {
            uint64_t y[22];
            int i;
            for (i = 0; i < 6; i++) y[i] = cur[5 * i + j];
            for (i = 6; i < 22; i++) y[i] = y[n1[i - 6]] ^ ROTL64(y[n2[i - 6]], m[i - 6]);
            cur[5 * 0 + j] = ROTL64(y[16], 0);
            cur[5 * 1 + j] = ROTL64(y[17], 15);
            cur[5 * 2 + j] = ROTL64(y[18], 48);
            cur[5 * 3 + j] = ROTL64(y[19], 63);
            cur[5 * 4 + j] = ROTL64(y[20], 6);
            cur[5 * 5 + j] = ROTL64(y[21], 33);
        }
        int i;
        for (i = 0; i < 6; i++) {
            uint64_t a0 = cur[5 * i + 0];
            uint64_t a1 = cur[5 * i + 1];
            uint64_t a2 = cur[5 * i + 2];
            uint64_t a3 = cur[5 * i + 3];
            uint64_t a4 = cur[5 * i + 4];
            cur[5 * i + 0] = a1 ^ (a0 & a1) ^ a2 ^ (a1 & a2) ^ a3 ^ (a3 & a4);
            cur[5 * i + 1] = ~0ULL ^ a1 ^ (a0 & a3) ^ (a1 & a3) ^ a4 ^ (a2 & a4);
            cur[5 * i + 2] = (a1 & a2) ^ a3 ^ (a2 & a3) ^ a4 ^ (a0 & a4);
            cur[5 * i + 3] = (a0 & a2) ^ (a1 & a3) ^ (a2 & a3) ^ a4;
            cur[5 * i + 4] = a0 ^ (a2 & a3) ^ (a1 & a4);
        }
        for (i = 0; i < 6; i++) {
            for (j = 0; j < 5; j++) {
                int n = 5 * i + j;
                int np = (10 * (n + 1) % 31) - 1;
                int ip = np / 5;
                int jp = np % 5;
                next[5 * ip + jp] = cur[n];
            }
        }
        for (int j2 = 0; j2 < 5; j2++) {
            int q = (b - 1) * 60 + 5 * t + j2;
            next[j2] ^= RC[q];
        }
        uint64_t *swap = cur;
        cur = next;
        next = swap;
    }
    if (cur != A) memcpy(A, tmp, 30 * sizeof(uint64_t));
}
static void TaiChi_P1(uint64_t state[30]) { P1920(state, 1); }
static void TaiChi_P2(uint64_t state[30]) { P1920(state, 2); }

static unsigned long long TaiChi_Padding(unsigned char *out, const unsigned char *msg, unsigned long long msg_len_bits) {
    unsigned long long bytes = msg_len_bits / 8;
    memcpy(out, msg, bytes);
    unsigned long long remain = msg_len_bits % 8;
    if (remain) out[bytes] = msg[bytes] & (0xff << (8 - remain));
    unsigned long long total_bits = msg_len_bits;
    out[total_bits / 8] |= (1 << (7 - total_bits % 8));
    total_bits++;
    while ((total_bits + 1) % RATE_BITS != 0) total_bits++;
    out[total_bits / 8] |= (1 << (7 - total_bits % 8));
    total_bits++;
    return (total_bits + 7) / 8;
}
static void TaiChi_Step(uint64_t S[RATE_WORDS], uint64_t L[CAPACITY_WORDS], uint64_t R[CAPACITY_WORDS], uint64_t M[RATE_WORDS], DOMAIN D) {
    uint64_t state[30] = {0};
    uint64_t A[RATE_WORDS] = {0};
    uint64_t B[CAPACITY_WORDS] = {0};
    int i;
    for (i = 0; i < RATE_WORDS; i++) state[i] = S[i] ^ M[i];
    for (i = 0; i < CAPACITY_WORDS; i++) state[RATE_WORDS + i] = L[i];
    if (D == DOMAIN_AB) state[29] ^= 0x01;
    else if (D == DOMAIN_FIN) state[29] ^= 0x02;
    else state[29] ^= 0x03;
    TaiChi_P1(state);
    memcpy(A, state, RATE_WORDS * sizeof(uint64_t));
    memcpy(B, state + RATE_WORDS, CAPACITY_WORDS * sizeof(uint64_t));
    for (i = 0; i < RATE_WORDS; i++) state[i] = A[i];
    for (i = 0; i < CAPACITY_WORDS; i++) state[RATE_WORDS + i] = R[i] ^ B[i];
    TaiChi_P2(state);
    for (i = 0; i < RATE_WORDS; i++) S[i] = state[i];
    for (i = 0; i < CAPACITY_WORDS; i++) { L[i] = state[RATE_WORDS + i]; R[i] = B[i] ^ state[RATE_WORDS + i]; }
}
static void TaiChi_Hash(const unsigned char *msg, unsigned long long msg_len_bits, unsigned char *digest) {
    uint64_t S[RATE_WORDS] = {0}, L[CAPACITY_WORDS] = {0}, R[CAPACITY_WORDS] = {0};
    unsigned long long padded_len = (msg_len_bits + 2 + RATE_BITS - 1) / RATE_BITS * RATE_BYTES;
    unsigned char *buf = (unsigned char *)malloc(padded_len);
    if (!buf) return;
    memset(buf, 0, padded_len);
    unsigned long long plen = TaiChi_Padding(buf, msg, msg_len_bits);
    uint64_t M[RATE_WORDS] = {0};
    unsigned long long blocks = plen / RATE_BYTES;
    for (unsigned long long i = 0; i < blocks; i++) {
        memcpy(M, buf + i * RATE_BYTES, RATE_BYTES);
        if (i == blocks - 1) TaiChi_Step(S, L, R, M, DOMAIN_FIN);
        else TaiChi_Step(S, L, R, M, DOMAIN_AB);
    }
    memcpy(digest, S, DIGEST_BYTES);
    free(buf);
}

/* ---- driver ---- */
static void hex(const unsigned char *b, int n) { for (int i = 0; i < n; i++) printf("%02x", b[i]); }

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "perm") == 0) {
        // Apply P1920(state, b) to a state seeded by splitmix64(seed); print input and output.
        uint64_t seed = strtoull(argv[2], 0, 0);
        int b = atoi(argv[3]);
        uint64_t st[30];
        uint64_t z = seed;
        for (int i = 0; i < 30; i++) {
            z += 0x9e3779b97f4a7c15ULL;
            uint64_t x = z;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            st[i] = x ^ (x >> 31);
        }
        printf("in :");
        for (int i = 0; i < 30; i++) printf(" %016llx", (unsigned long long)st[i]);
        printf("\n");
        P1920(st, b);
        printf("out:");
        for (int i = 0; i < 30; i++) printf(" %016llx", (unsigned long long)st[i]);
        printf("\n");
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "vectors") == 0) {
        // For seed=0..N-1 and b=1,2: seed a state via splitmix64, apply P1920, dump output words.
        int N = argc >= 3 ? atoi(argv[2]) : 200;
        for (int s = 0; s < N; s++) {
            for (int b = 1; b <= 2; b++) {
                uint64_t st[30];
                uint64_t z = (uint64_t)s;
                for (int i = 0; i < 30; i++) {
                    z += 0x9e3779b97f4a7c15ULL;
                    uint64_t x = z;
                    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
                    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
                    st[i] = x ^ (x >> 31);
                }
                P1920(st, b);
                printf("%d %d", s, b);
                for (int i = 0; i < 30; i++) printf(" %016llx", (unsigned long long)st[i]);
                printf("\n");
            }
        }
        return 0;
    }
    // Default: hash a few messages and print digests.
    const char *tests[] = {"", "abc", "The quick brown fox jumps over the lazy dog"};
    for (int t = 0; t < 3; t++) {
        unsigned char dig[DIGEST_BYTES];
        unsigned long long bitlen = (unsigned long long)strlen(tests[t]) * 8;
        TaiChi_Hash((const unsigned char *)tests[t], bitlen, dig);
        printf("msg=\"%s\" (%llu bits)\n  digest=", tests[t], bitlen);
        hex(dig, DIGEST_BYTES);
        printf("\n");
    }
    // Also a 200-byte incremental message, for a multi-block check.
    unsigned char big[200];
    for (int i = 0; i < 200; i++) big[i] = (unsigned char)(i * 37 + 11);
    unsigned char dig[DIGEST_BYTES];
    TaiChi_Hash(big, 200 * 8, dig);
    printf("msg=200-byte-ramp (1600 bits)\n  digest=");
    hex(dig, DIGEST_BYTES);
    printf("\n");
    return 0;
}
