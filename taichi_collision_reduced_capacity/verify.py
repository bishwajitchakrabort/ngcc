#!/usr/bin/env python3
"""Independent check of a reduced-capacity TaiChi collision.

Reimplements the real P1920 permutation (forward only) and the reduced IFS mode from
scratch, then hashes the two raw messages written by `taichi_real ... --out PREFIX`
(PREFIX_M.hex, PREFIX_Mprime.hex, each "bitlen hex") and checks the digests collide.

usage: verify.py c PREFIX_M.hex PREFIX_Mprime.hex
       (c is the capacity used in the run, e.g. 40)
"""
import sys

MASK = (1 << 64) - 1
def rotl(x, n): return ((x << n) | (x >> (64 - n))) & MASK if n else x

N1 = [0, 1, 2, 4, 3, 5, 6, 7, 13, 14, 15, 16, 15, 18, 11, 20]
N2 = [1, 2, 4, 3, 5, 6, 0, 11, 8, 9, 10, 12, 17, 14, 19, 13]
MM = [8, 0, 0, 0, 0, 0, 23, 62, 35, 14, 48, 1, 57, 63, 58, 22]
OROT = [0, 15, 48, 63, 6, 33]

def make_RC():
    p = [0, 2, 3, 6, 13, 28, 59]
    s = [1, 0, 0, 0, 0, 0, 0]
    RC = []
    for _ in range(120):
        C = 0
        for i in range(7):
            if s[i]: C ^= (1 << p[i])
        RC.append(C)
        n = [s[6], s[0] ^ s[6], s[1], s[2], s[3], s[4], s[5]]
        s = n
    return RC
RC = make_RC()
K2 = [(10 * (n + 1) % 31) - 1 for n in range(30)]

def P1920(A, b):
    cur = list(A)
    for t in range(12):
        nxt = [0] * 30
        for j in range(5):
            y = [cur[5 * i + j] for i in range(6)] + [0] * 16
            for i in range(6, 22):
                y[i] = y[N1[i - 6]] ^ rotl(y[N2[i - 6]], MM[i - 6])
            for i in range(6):
                cur[5 * i + j] = rotl(y[16 + i], OROT[i])
        for i in range(6):
            a0, a1, a2, a3, a4 = (cur[5 * i + k] for k in range(5))
            cur[5 * i + 0] = a1 ^ (a0 & a1) ^ a2 ^ (a1 & a2) ^ a3 ^ (a3 & a4)
            cur[5 * i + 1] = MASK ^ a1 ^ (a0 & a3) ^ (a1 & a3) ^ a4 ^ (a2 & a4)
            cur[5 * i + 2] = (a1 & a2) ^ a3 ^ (a2 & a3) ^ a4 ^ (a0 & a4)
            cur[5 * i + 3] = (a0 & a2) ^ (a1 & a3) ^ (a2 & a3) ^ a4
            cur[5 * i + 4] = a0 ^ (a2 & a3) ^ (a1 & a4)
        for n in range(30):
            nxt[K2[n]] = cur[n]
        for j in range(5):
            nxt[j] ^= RC[(b - 1) * 60 + 5 * t + j]
        cur = nxt
    return cur

def taichi_reduced(msg, bits, c):
    R = 1920 - c
    RB = R // 8
    capshift = 64 - c
    lowmask = (1 << capshift) - 1 if capshift else 0
    cmask = (1 << c) - 1
    buf = bytearray(msg)
    # 10*1 padding, RATE_BITS = R
    total = bits
    buf += bytes((len(buf) * 8 + 8 - total) // 8 + 64)  # ensure room
    buf = bytearray(msg) + bytearray((max(0, (((bits + 2 + R - 1) // R) * RB) - len(msg)) + 8))
    if bits % 8:
        buf[bits // 8] &= (0xff << (8 - bits % 8)) & 0xff
    buf[total // 8] |= 1 << (7 - total % 8); total += 1
    while (total + 1) % R != 0: total += 1
    buf[total // 8] |= 1 << (7 - total % 8); total += 1
    plen = (total + 7) // 8
    nb = plen // RB
    def cap_get(w): return (w >> capshift) if capshift else w
    def cap_set(w, v): return ((w & lowmask) | ((v & cmask) << capshift)) if capshift else (v & cmask)
    S = [0] * 30; L = 0; Rr = 0
    for i in range(nb):
        blk = buf[i * RB:(i + 1) * RB]
        M = [0] * 30
        for bb in range(len(blk)):
            M[bb // 8] |= blk[bb] << (8 * (bb % 8))
        tau = 2 if i == nb - 1 else 1
        st = [S[k] ^ M[k] for k in range(30)]
        st[29] = cap_set(st[29], L ^ tau)
        st = P1920(st, 1)
        Bc = cap_get(st[29])
        st[29] = cap_set(st[29], Rr ^ Bc)
        st = P1920(st, 2)
        C = cap_get(st[29])
        S = list(st); S[29] &= lowmask
        L = C; Rr = (Bc ^ C) & cmask
    dig = bytearray()
    for k in range(c // 8):
        dig += bytes([(S[k // 8] >> (8 * (k % 8))) & 0xff])
    return S, L, Rr, dig

def load(fn):
    bits, hx = open(fn).read().split()
    return int(bits), bytes.fromhex(hx)

def main():
    c = int(sys.argv[1])
    b1, m1 = load(sys.argv[2])
    b2, m2 = load(sys.argv[3])
    S1, L1, R1, d1 = taichi_reduced(m1, b1, c)
    S2, L2, R2, d2 = taichi_reduced(m2, b2, c)
    print(f"c = {c}")
    print(f"M  : {b1} bits, digest {d1.hex()}")
    print(f"M' : {b2} bits, digest {d2.hex()}")
    print(f"final S equal : {S1 == S2}")
    print(f"final L equal : {L1 == L2}")
    print(f"R != R' (expected): {R1 != R2}")
    ok = (m1 != m2) and (S1 == S2) and (d1 == d2)
    print("COLLISION VERIFIED" if ok else "NOT A COLLISION")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
