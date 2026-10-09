#!/usr/bin/env python3
# ============================================================================
# chash1024_2ndpreimage_toy.py
#
# Toy-scale, self-contained demonstration of the backward-tree second-preimage
# attack on the Rocket-JH / CTR-Func mode of C Hash (C-Hash-1024).
#
# WHAT THIS IS
#   A scaled-down model of the mode. The 1536-bit C-Engine permutation is
#   replaced by a small balanced Feistel permutation on w = b + t bits, the
#   effective state b is 16..24 instead of 1472, and the counter t is a few
#   bits. The mode logic (CTR-Func feed-forward, JH message injection,
#   per-position counters, FIL finalization) is identical to the spec. The
#   purpose is to confirm that the attack and its complexity formula, derived
#   for the real parameters, behave as predicted when shrunk to a size that
#   runs on a laptop.
#
# WHAT IT DEMONSTRATES
#   1. demo  : takes ONE challenge message (one hash chain), finds a DIFFERENT
#              message with the same digest, and verifies it. This is a genuine
#              second preimage.
#   2. growth: the backward tree grows as B_j = (N / 2^(b/2))^j, for several b.
#   3. work  : the cost per second preimage follows
#                 (d+1) * 2^((d*(b/2) + b)/(d+1))
#              and decreases as the tree depth d increases.
#
# HOW TO RUN
#   python3 chash1024_2ndpreimage_toy.py            # demo + quick self-checks
#   python3 chash1024_2ndpreimage_toy.py --demo     # just the worked example
#   python3 chash1024_2ndpreimage_toy.py --growth   # just the growth law
#   python3 chash1024_2ndpreimage_toy.py --work     # just the work scaling
#   python3 chash1024_2ndpreimage_toy.py --b 24 --d 3 --demo
#
# MAPPING TO THE REAL INSTANCE
#   real C-Hash-1024:  b (state) = 1472, t (counter) = 64, perm width = 1536,
#   message block = b/2 = 736, 2nd-preimage claim = 2^1024.
#   The same formula gives 2^983 at d=2 and 2^922 at d=3, both below 2^1024.
#
# No third-party dependencies. Python 3.8+.
# ============================================================================
import argparse
import math
import random

MASK64 = (1 << 64) - 1


def splitmix(z):
    """A 64-bit mixing function, used only to build the toy round function."""
    z = (z + 0x9E3779B97F4A7C15) & MASK64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return z ^ (z >> 31)


class Perm:
    """Balanced Feistel permutation on w bits, standing in for the C-Engine P.
    Invertible, so the same object can model both CTR-Perm and CTR-Func."""

    def __init__(self, w, seed, rounds=12):
        assert w % 2 == 0
        self.w = w
        self.half = w // 2
        self.hmask = (1 << self.half) - 1
        rng = random.Random(seed)
        self.keys = [rng.getrandbits(64) for _ in range(rounds)]

    def fwd(self, v):
        L = v >> self.half
        R = v & self.hmask
        for k in self.keys:
            L, R = R, L ^ (splitmix(R ^ k) & self.hmask)
        return (L << self.half) | R


class ToyHash:
    """Rocket-JH in CTR-Func mode on a (b+t)-bit permutation.

    State b bits, two halves of b/2. Message block b/2 bits.
    VIL round (matches Alg.1 with P'_ctr(x) = P_ctr(x) XOR x):
        x_i   = h_{i-1} XOR (m_i || 0)          inject message into left half
        G(x)  = LMB(b, P(x || ctr_i)) XOR x     CTR-Func feed-forward, G = F+id
        h_i   = G(x_i) XOR (0 || m_i)           inject message into right half
    FIL finalization (one block, as in C-Hash-1024, z=1):
        digest = ( LMB(b, P(h_l || ctr_fil)) XOR h_l ) truncated to hashlen
    """

    def __init__(self, b, t, seed=1, hashlen=None):
        assert b % 2 == 0 and t >= 2
        self.b, self.t = b, t
        self.w = b + t
        self.hb = b // 2
        self.bmask = (1 << b) - 1
        self.hmask = (1 << self.hb) - 1
        self.tmask = (1 << t) - 1
        self.P = Perm(self.w, seed)
        self.hashlen = hashlen if hashlen is not None else b
        self.ctr_fil = (1 << (t - 1)) | 1          # FIL lane, distinct from VIL
        self.h0 = splitmix(seed ^ 0xABCDEF) & self.bmask

    def Ftrunc(self, x, ctr):
        """LMB(b, P(x || ctr)): top b bits of the permutation output."""
        y = self.P.fwd(((x & self.bmask) << self.t) | (ctr & self.tmask))
        return y >> self.t

    def ctr_vil(self, pos):
        """Per-position VIL counter (lane bit 0). Distinct per block index."""
        return pos & self.tmask

    def rnd(self, h_prev, m, ctr):
        """One CTR-Func VIL round."""
        x = h_prev ^ ((m & self.hmask) << self.hb)
        G = self.Ftrunc(x, ctr) ^ x
        return (G ^ (m & self.hmask)) & self.bmask

    def chain(self, blocks):
        """Absorb a list of b/2-bit blocks, return the final VIL state."""
        h = self.h0
        for i, m in enumerate(blocks, start=1):
            h = self.rnd(h, m, self.ctr_vil(i))
        return h

    def fil(self, h_last):
        d = (self.Ftrunc(h_last, self.ctr_fil) ^ h_last) & self.bmask
        return d & ((1 << self.hashlen) - 1)

    def digest(self, blocks):
        return self.fil(self.chain(blocks))

    def backstep(self, target, x, ctr):
        """Backward step. Given a target state and a free input x, return a
        predecessor (h_prev, m) with rnd(h_prev, m, ctr) == target, or None.
        Valid iff G(x) agrees with target on the left (high) half."""
        G = self.Ftrunc(x, ctr) ^ x
        if (G >> self.hb) != (target >> self.hb):
            return None
        m = (G & self.hmask) ^ (target & self.hmask)
        h_prev = x ^ ((m & self.hmask) << self.hb)
        return h_prev, m


# ---------------------------------------------------------------------------
# The attack
# ---------------------------------------------------------------------------
def forward_states(H, depth, F, rng):
    """~F reachable states at the given depth, with the message prefix that
    reaches each. Varies the last two blocks so the diversity ceiling is 2^b
    rather than 2^(b/2) (a single block is only b/2 bits)."""
    k = 2 if depth >= 2 else 1
    stem = [rng.getrandbits(H.hb) for _ in range(depth - k)]
    h_stem = H.chain(stem) if stem else H.h0
    out = {}
    guard = 0
    while len(out) < F and guard < 3 * F + 16:
        guard += 1
        varied = [rng.getrandbits(H.hb) for _ in range(k)]
        h = h_stem
        for j in range(k):
            h = H.rnd(h, varied[j], H.ctr_vil(depth - k + 1 + j))
        out[h] = stem + varied
    return out


def attack_once(H, L, d, F, N, rng):
    """One attempt on one challenge chain of length L, tree depth d.
    Returns (challenge_blocks, forged_blocks) for the first verified second
    preimage, or None."""
    chal = [rng.getrandbits(H.hb) for _ in range(L)]
    Y = H.chain(chal)                       # target = final VIL state of challenge
    chal_dig = H.fil(Y)

    # Backward tree from Y over positions L, L-1, ..., L-d+1.
    # Each node carries the tail of blocks that drives it to Y.
    level = {Y: []}
    for j in range(d):
        ctr = H.ctr_vil(L - j)
        tbl = {}
        for st, tail in level.items():
            tbl.setdefault(st >> H.hb, []).append((st, tail))
        nxt = {}
        for _ in range(N):
            x = rng.getrandbits(H.b)
            G = H.Ftrunc(x, ctr) ^ x
            bucket = tbl.get(G >> H.hb)
            if bucket:
                for st, tail in bucket:
                    m = (G & H.hmask) ^ (st & H.hmask)
                    nxt[x ^ ((m & H.hmask) << H.hb)] = [m] + tail
        level = nxt
        if not level:
            return None
    leaves = level                          # states at depth L-d, with d-block tails

    # Forward side: F reachable states at depth L-d. Match on the full b bits.
    fwd = forward_states(H, L - d, F, rng)
    for st, tail in leaves.items():
        if st in fwd:
            forged = fwd[st] + tail
            if forged != chal and len(forged) == L and H.digest(forged) == chal_dig:
                return chal, forged
    return None


def hexblocks(H, blocks):
    width = (H.hb + 3) // 4
    return " ".join(f"{m & H.hmask:0{width}x}" for m in blocks)


def demo(b=20, t=8, d=2, L=8, seed=2026, bump=3.0):
    """One challenge chain -> one verified second preimage."""
    H = ToyHash(b, t, seed=seed)
    s = b // 2
    logT = (d * s + b) / (d + 1)
    T = max(2, int(bump * 2 ** logT))
    F = N = T
    print(f"[demo] toy CTR-Func  b(state)={b}  t(ctr)={t}  perm={b+t}  "
          f"block={s}  depth d={d}")
    print(f"[demo] T = 2^{logT:.2f} (x{bump} margin),  challenge length L={L}")
    rng = random.Random(seed)
    for attempt in range(1, 61):
        r = attack_once(H, L, d, F, N, rng)
        if r:
            chal, forged = r
            dc, df = H.digest(chal), H.digest(forged)
            print(f"[demo] found on attempt {attempt}")
            print(f"       challenge M  = {hexblocks(H, chal)}")
            print(f"       forged    M' = {hexblocks(H, forged)}")
            print(f"       H(M)  = {dc:0{(H.hashlen+3)//4}x}")
            print(f"       H(M') = {df:0{(H.hashlen+3)//4}x}")
            print(f"       M != M' : {chal != forged};   "
                  f"same length : {len(chal) == len(forged)};   "
                  f"same digest : {dc == df}")
            print("[demo] VERIFIED second preimage" if (chal != forged and dc == df)
                  else "[demo] FAILED")
            return
    print("[demo] no second preimage within the attempt budget (raise --bump)")


def growth(bs=(16, 20, 24, 28), dmax=6, rho=4, trials=10):
    print(f"[growth] B_j vs (N/2^(b/2))^j   N = {rho}*2^(b/2)   {trials} trials")
    print("   b | " + " ".join(f"B{j+1:<7d}" for j in range(dmax)))
    print("     | " + " ".join(f"{rho**(j+1):<8d}" for j in range(dmax)) + " (predicted)")
    for b in bs:
        H = ToyHash(b, 8, seed=1000 + b)
        s = b // 2
        N = rho * (1 << s)
        acc = [0.0] * dmax
        for tr in range(trials):
            rng = random.Random((b << 10) ^ tr ^ 0x5151)
            level = {rng.getrandbits(b): None}
            for j in range(dmax):
                ctr = (200 + j) & H.tmask
                tbl = {}
                for st in level:
                    tbl.setdefault(st >> H.hb, []).append(st)
                nxt = {}
                for _ in range(N):
                    x = rng.getrandbits(b)
                    G = H.Ftrunc(x, ctr) ^ x
                    bucket = tbl.get(G >> H.hb)
                    if bucket:
                        for st in bucket:
                            m = (G & H.hmask) ^ (st & H.hmask)
                            nxt[x ^ ((m & H.hmask) << H.hb)] = None
                acc[j] += len(nxt)
                level = nxt
                if not level:
                    break
        print(f" {b:3d} | " + " ".join(f"{acc[j]/trials:<8.1f}" for j in range(dmax)))


def work(b=24, ds=(1, 2, 3, 4), trials=16, L=8):
    print(f"[work] cost per 2nd-preimage vs (d+1)*2^((d*b/2+b)/(d+1))   b={b}")
    print("   d |  T=2^logT | verified/trials | work obs | work pred")
    for d in ds:
        H = ToyHash(b, 8, seed=500 + d)
        s = b // 2
        logT = (d * s + b) / (d + 1)
        T = int(round(2 ** logT))
        F = N = T
        rng = random.Random(77 + d)
        ok = 0
        for _ in range(trials):
            if attack_once(H, L, d, F, N, rng):
                ok += 1
        pred = logT + math.log2(d + 1)
        obs = math.log2((d + 1) * T * trials / ok) if ok else float("nan")
        print(f"   {d} | 2^{logT:5.2f}   | {ok:3d}/{trials:<11d} | "
              f"2^{obs:6.2f} | 2^{pred:.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Toy C-Hash-1024 second-preimage attack")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--growth", action="store_true")
    ap.add_argument("--work", action="store_true")
    ap.add_argument("--b", type=int, default=20)
    ap.add_argument("--t", type=int, default=8)
    ap.add_argument("--d", type=int, default=2)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--bump", type=float, default=3.0)
    args = ap.parse_args()

    run_all = not (args.demo or args.growth or args.work)
    if args.demo or run_all:
        demo(b=args.b, t=args.t, d=args.d, seed=args.seed, bump=args.bump)
        print()
    if args.growth or run_all:
        growth()
        print()
    if args.work or run_all:
        work()
