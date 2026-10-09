"""Toy check of the FChain palindrome prefix attack (WChain mode, ideal-ish P_B).

Mode as in the WChain spec (Alg. 1, eqs. 6-9): H_{-1}=0, H_0=IV,
H_i = CF(H_{i-1}, M_i) xor H_{i-2},  H_{l+1} = CF(H_l, Sigma),  Z = Trunc_h.
CF(V,B) = P_B(V) xor V, with P_B a keyed permutation (4-round Feistel here).
10* padding, no length field, Sigma = xor of padded blocks.
"""
import hashlib
import math
import random

B_BITS = 20            # toy state b
HALF = B_BITS // 2
MASK = (1 << B_BITS) - 1
HMASK = (1 << HALF) - 1
MU = 2 * B_BITS        # block length, mu = 2b as in WChain
H_BITS = 16            # truncated digest
ROUNDS = 4

calls = {"fwd": 0, "inv": 0}


def rf(key, r, x):
    d = hashlib.blake2b(key.to_bytes(8, "big") + bytes([r]) + x.to_bytes(2, "big"),
                        digest_size=4).digest()
    return int.from_bytes(d, "big") & HMASK


def P(key, v):
    calls["fwd"] += 1
    l, r = v >> HALF, v & HMASK
    for i in range(ROUNDS):
        l, r = r, l ^ rf(key, i, r)
    return (l << HALF) | r


def Pinv(key, v):
    calls["inv"] += 1
    l, r = v >> HALF, v & HMASK
    for i in reversed(range(ROUNDS)):
        l, r = r ^ rf(key, i, l), l
    return (l << HALF) | r


def CF(v, blk):
    return P(blk, v) ^ v


def pad(bits):
    k = (-(len(bits) + 1)) % MU
    s = bits + "1" + "0" * k
    return [int(s[i:i + MU], 2) for i in range(0, len(s), MU)]


def wchain_mode(bits, iv):
    blocks = pad(bits)
    sigma = 0
    for m in blocks:
        sigma ^= m
    hm2, hm1 = 0, iv
    for m in blocks:
        hm2, hm1 = hm1, CF(hm1, m) ^ hm2
    return CF(hm1, sigma) >> (B_BITS - H_BITS)


def birthday(left, right, n):
    """Find (X, Y) with left(X) == right(Y); n candidates per side."""
    table = {}
    for _ in range(n):
        x = random.getrandbits(MU)
        table[left(x)] = x
    for _ in range(n):
        y = random.getrandbits(MU)
        t = right(y)
        if t in table:
            return table[t], y
    return None


def blk_bits(x):
    return format(x, "0%db" % MU)


def check(prefix_bits, iv, trials=2000):
    for _ in range(trials):
        m = "".join(random.choice("01") for _ in range(random.randrange(0, 5 * MU)))
        if wchain_mode(prefix_bits + m, iv) != wchain_mode(m, iv):
            return False
    return True


random.seed(1)
n = 1 << (B_BITS // 2 + 1)

# Case 1: IV = 0, single palindrome XYYX with P_X(0) = P_Y^{-1}(0).
calls.update(fwd=0, inv=0)
X, Y = None, None
while X is None:
    r = birthday(lambda x: P(x, 0), lambda y: Pinv(y, 0), n)
    if r:
        X, Y = r
cost1 = calls["fwd"] + calls["inv"]
pref1 = "".join(blk_bits(b) for b in (X, Y, Y, X))
print("IV=0  XYYX     : search calls = 2^%.1f (b/2 = %d), Hash(XYYX||M)==Hash(M) on all trials: %s"
      % (math.log2(cost1), B_BITS // 2, check(pref1, 0)))

# Case 2: IV != 0. One palindrome only swaps (0,IV) -> (IV,0); two palindromes close the loop.
IV = random.getrandbits(B_BITS) | 1
print("IV!=0 XYYX only: Hash equal on all trials: %s (expected False)" % check(pref1, IV, trials=200))
calls.update(fwd=0, inv=0)
r1 = None
while r1 is None:
    r1 = birthday(lambda x: P(x, IV) ^ IV, lambda y: Pinv(y, IV), n)
r2 = None
while r2 is None:
    r2 = birthday(lambda z: P(z, 0) ^ IV, lambda w: Pinv(w, 0), n)
cost2 = calls["fwd"] + calls["inv"]
(X1, Y1), (Z1, W1) = r1, r2
pref2 = "".join(blk_bits(b) for b in (X1, Y1, Y1, X1, Z1, W1, W1, Z1))
print("IV!=0 XYYX ZWWZ: search calls = %d (vs %d for one match), Hash equal on all trials: %s"
      % (cost2, cost1, check(pref2, IV)))
