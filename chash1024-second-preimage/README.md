# C-Hash-1024 second-preimage attack

A multi-level backward-tree second-preimage attack on **C-Hash-1024**, the
Rocket-JH construction in CTR-Func mode on the 1536-bit C-Engine permutation
(C Hash design document, 30 June 2026).

The attack finds a second preimage in about `2^983` permutation calls with
`2^491` memory (tree depth `d=2`), down to about `2^746` at the memory-heavy
end. All sit below the claimed `2^1024` second-preimage resistance (Table 6,
Appendix C.7.2). The cost is governed by the effective state width through the
`n/2` message-injection rate and does not depend on the counter width `w`.
Preimage resistance, collision resistance, and every C-Hash-512 claim are
unaffected.

This folder includes a toy-scale, self-contained verification. It shrinks the
state to 16-24 bits and replaces the C-Engine with a small Feistel permutation,
keeping the mode logic identical, so the attack and its complexity formula can
be checked on a laptop in a few minutes.

## Contents

- `README.md` - this file, including the full analysis
- `chash1024_2ndpreimage_toy.py` - the verification script (no dependencies)

## Requirements

Python 3.8 or later. Standard library only.

## Run

```
python3 chash1024_2ndpreimage_toy.py            # demo + growth + work checks
python3 chash1024_2ndpreimage_toy.py --demo     # one-chain worked example only
python3 chash1024_2ndpreimage_toy.py --growth   # tree-growth law only
python3 chash1024_2ndpreimage_toy.py --work     # work scaling only
python3 chash1024_2ndpreimage_toy.py --b 24 --d 3 --demo
```

In the code, `b` is the effective state `n` of the analysis below and `t` is the
counter `w`.

---

## Notation

We follow the spec. Effective state `n = 1472`, counter `w = 64`, permutation
width `b = n + w = 1536`, message block `n/2 = 736`. A state is written as two
halves `(left || right)`, each `n/2` bits. The retained map is
`P_ctr(x) = LMB(n, P(x || ctr))`, and the CTR-Func round applies the feed-forward
`P'_ctr(x) = P_ctr(x) ⊕ x`. The VIL phase uses a per-position counter `ctr_i`;
the FIL phase uses a fixed, lane-separated counter `ctr_fil`. A challenge message
is `M'` with `l` blocks and digest `y`.

## The round map collapses to one call plus two injections

From Algorithm 1 with the CTR-Func definition, one VIL round is

```
x_i = h_{i-1} ⊕ (m_i || 0)
h_i = P'_ctr(x_i) ⊕ (0 || m_i) = F_i(x_i) ⊕ x_i ⊕ (0 || m_i)
```

with `F_i(x) = LMB(n, P(x || ctr_i))`. Writing `G_i = F_i ⊕ id`, this is

```
h_i = G_i(x_i) ⊕ (0 || m_i).
```

So the left half of `h_i` equals `G_i(x_i)` on its left half, an `n/2`-bit value,
and the right half equals `G_i(x_i)_R ⊕ m_i`. The message in the right half is
free.

## The backward step is a forward-only n/2-bit condition

Fix a target state `h` and a counter `ctr`. Pick any input `x` and compute
`G(x) = F(x) ⊕ x`. Then `x` is a predecessor-generating input for `h` exactly
when

```
G(x)_L = h_L          (an n/2-bit condition)
```

in which case `m = G(x)_R ⊕ h_R` and the predecessor is `h_prev = x ⊕ (m || 0)`.
One verifies `rnd(h_prev, m, ctr) = h` by construction. This uses forward calls
to `P` only. It does not need `P_ctr` to be invertible, so the CTR-Func
feed-forward, whose whole purpose is to block backward queries, does not block
this step. A fresh `x` meets the condition for a given target with probability
`2^-(n/2)`, and the same `G(x)` can be tested against many targets at once by
indexing them on their left half.

## The backward tree

Grow a tree of predecessors backward from the final VIL state `Y = h_l` of the
challenge, over the counters `ctr_l, ctr_{l-1}, ..., ctr_{l-d+1}`.

```
Level 1 (ctr_l):       N forward probes, keep predecessors of Y.
                       B_1 = N / 2^(n/2).
Level j (ctr_{l-j+1}): N forward probes, each matched on the left half against
                       all B_{j-1} current nodes.  B_j = N * B_{j-1} / 2^(n/2).
After d levels:        B_d = N^d / 2^(d*n/2) leaves, each a state at depth l-d
                       carrying a d-block tail that drives it to Y.
Forward side:          F reachable states at depth l-d, from the IV. Match the
                       leaves against them on the full n bits: F * B_d / 2^n
                       expected collisions.
```

A collision gives a message of exactly `l` blocks: an `(l-d)`-block prefix from
the IV to the matched state, then the `d`-block tail to `Y`, then nothing further
because `Y` is already the final VIL state. Appending the challenge's own FIL
gives the same digest `y`, since FIL reads only `h_l` and its fixed counter. The
length is unchanged, so the padding is identical, and the message differs from
the challenge. This is a genuine second preimage.

## Complexity

With `N = F = T` and success `F * B_d / 2^n = 1`,

```
log2 T      = (d*(n/2) + n)/(d+1) = (n/2)*(d+2)/(d+1)
log2 time   = log2 T + log2(d+1)
log2 memory = (n/2) * d/(d+1)
floor (d -> infinity) = n/2
```

For `n = 1472`:

| d | log2 time | log2 memory |
|---|-----------|-------------|
| 1 | 1105 | 368 |
| 2 | **983** | 491 |
| 3 | 922 | 552 |
| 4 | 886 | 589 |
| 8 | 821 | 654 |
| ~509 (opt) | 746 | 735 |

`d=1` is the ordinary meet-in-the-middle and stays above the claim. Every
`d >= 2` is below `2^1024`. The floor is `2^(n/2) = 2^736`, reached only with
about `2^735` memory. The practical point is `d=2`: `2^983` time, `2^491`
memory.

The floor equals `2^(n-r)` for an injection rate `r`, and here `r = n/2`. The
736-bit absorption rate the design advertises is exactly what pins the floor at
`2^736`.

## Scope

The attack is specific to second preimage on C-Hash-1024.

- **C-Hash-512** (CTR-Perm, claim 512): the tree costs at least `2^746`; the
  cheaper CTR-Perm route, an inverse meet-in-the-middle at `2^((n+w)/2) = 2^768`,
  also stays far above 512. Not threatened.
- **Collisions** (claims 256 and 512): the `2^768` inverse route is the cheapest
  internal collision on the shared permutation. Both claims are below it. Not
  threatened.
- **Preimage** (claims 512 and 1024): unaffected. With no target state in hand
  the adversary must invert the FIL, which is `2^hashlen`. Appendix C.7.1 holds.

Note that CTR-Func was chosen precisely to deny the `2^768` inverse route of
CTR-Perm. It does. The forward tree does not use that route and is not blocked
by it.

## Where the proof misses it

Two places, both the `d=1` truncation of the same tree.

Appendix C.5.1, the single-round damping argument, shows that `w` right-half
collisions cannot be advanced together, so the multiplicity dies in one round.
That is about advancing one query into many states. The tree never does that. It
matches one forward value `G(x)` against all current targets on the `n/2`-bit
left half, a hash-table lookup, so every query can land on any live target. The
correct count of states reachable at depth `d` is `(q / 2^(n/2))^d`, not the
`q / 2^(d*n/2)` implied by treating depth as independent per-query conditions.

Appendix C.7.2, Strategy 1, bounds the iterative-collision route by
`O(q^2 (l+l')) / 2^(3n/2)`. That balances at `q = 2^(3n/4) = 2^1104`, which is
exactly the `d=1` tree. The multi-level tree, `d >= 2`, is unaccounted for, and
it is what drops under `2^1024`.

## What the script checks

The code replaces the 1536-bit C-Engine with a small balanced Feistel
permutation on `n + w` bits and shrinks `n` to 16-24 bits, keeping the mode
logic identical: the CTR-Func feed-forward, the JH message injection, the
per-position VIL counters, and the one-call FIL.

1. **demo** takes one challenge message, builds its hash chain, and finds a
   different message of the same length with the same digest:

   ```
   [demo] found on attempt 1
          challenge M  = 079 147 202 3ce 20c 296 370 3cf
          forged    M' = 069 097 1d1 19d 1ff 2ec 39c 2dc
          H(M)  = a641a
          H(M') = a641a
          M != M' : True;   same length : True;   same digest : True
   [demo] VERIFIED second preimage
   ```

2. **growth** measures the tree size `B_j` against `(N / 2^(n/2))^j`. With
   `N = 4*2^(n/2)` the columns track `4, 16, 64, 256, 1024, 4096` across every
   `n`, confirming the exponent depends on `n` only through `n/2`, not on `w`.

3. **work** measures the cost per second preimage against
   `(d+1) * 2^((d*n/2 + n)/(d+1))`, showing it fall as the depth `d` rises.

One implementation note: the forward side must vary at least two message blocks
to reach `F` distinct states, since one block is only `n/2` bits and the forward
list would otherwise saturate at `2^(n/2)`. This costs nothing in the attack.

## AI disclosure

The cryptanalysis here (the attack idea, the round-map reduction, the
backward-tree construction, the complexity analysis, and this code) was produced
by an AI assistant (Claude, Anthropic). A human directed the work, reviewed the
output, and ran the code. The empirical results are at toy scale. The `n = 1472`
figures are closed-form, not a full-scale computation. Verify independently; the
script is provided for that.

## Author

Bishwajit Chakraborty, bishu.math.ynwa@gmail.com
