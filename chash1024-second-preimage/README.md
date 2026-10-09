# C-Hash-1024 second-preimage attack (toy verification)

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

This folder holds a toy-scale, self-contained verification. It shrinks the
state to 16-24 bits and replaces the C-Engine with a small Feistel permutation,
keeping the mode logic identical, so the attack and its complexity formula can
be checked on a laptop in a few minutes.

## Contents

- `chash1024_2ndpreimage_toy.py` - the verification script (no dependencies)
- `README.md` - this file

## Requirements

Python 3.8 or later. Standard library only.

## Run

```
python3 chash1024_2ndpreimage_toy.py            # demo + growth + work checks
python3 chash1024_2ndpreimage_toy.py --demo     # one-chain worked example only
python3 chash1024_2ndpreimage_toy.py --growth    # tree-growth law only
python3 chash1024_2ndpreimage_toy.py --work      # work scaling only
python3 chash1024_2ndpreimage_toy.py --b 24 --d 3 --demo
```

In the code, `b` is the effective state `n` of the write-up and `t` is the
counter `w`.

## What it checks

1. **demo** takes one challenge message, builds its hash chain, and finds a
   different message of the same length with the same digest. It prints both
   messages and their equal digests:

   ```
   [demo] found on attempt 1
          challenge M  = 079 147 202 3ce 20c 296 370 3cf
          forged    M' = 069 097 1d1 19d 1ff 2ec 39c 2dc
          H(M)  = a641a
          H(M') = a641a
          M != M' : True;   same length : True;   same digest : True
   [demo] VERIFIED second preimage
   ```

2. **growth** measures the backward-tree size `B_j` and compares it to
   `(N / 2^(n/2))^j`. With `N = 4*2^(n/2)` the columns track `4, 16, 64, 256,
   1024, 4096` across every `n`, confirming the exponent depends on `n` only
   through `n/2`, not on `w`.

3. **work** measures the cost per second preimage and compares it to
   `(d+1) * 2^((d*n/2 + n)/(d+1))`, showing the cost fall as the depth `d`
   rises.

## Complexity at the real parameters (n = 1472)

| d | log2 time | log2 memory |
|---|-----------|-------------|
| 1 | 1105 | 368 |
| 2 | **983** | 491 |
| 3 | 922 | 552 |
| 4 | 886 | 589 |
| 8 | 821 | 654 |
| ~509 (opt) | 746 | 735 |

`d=1` is the ordinary meet-in-the-middle and stays above the claim. Every
`d >= 2` is below `2^1024`. The floor is `2^(n/2) = 2^736`. The figures for
`n = 1472` are closed-form evaluations of the formula; the script verifies the
same formula at toy scale.

## Why the proof misses it

Appendix C.5.1 (single-round damping) bounds the states reachable at depth `d`
as if depth added independent per-query conditions. The tree matches each
forward value against all current targets on the `n/2`-bit left half, so every
query can land on any live target, and the correct count is `(q / 2^(n/2))^d`.
Appendix C.7.2 Strategy 1 bounds the iterative-collision route by
`O(q^2 (l+l')) / 2^(3n/2)`, which balances at `q = 2^(3n/4) = 2^1104`: exactly
the `d=1` case. The multi-level tree, `d >= 2`, is what drops under `2^1024`.

## AI disclosure

The cryptanalysis here (the attack idea, the round-map reduction, the
backward-tree construction, the complexity analysis, and this code) was
produced by an AI assistant (Claude, Anthropic). A human directed the work,
reviewed the output, and ran the code. The empirical results are at toy scale.
The `n = 1472` figures are closed-form, not a full-scale computation. Verify
independently; the script is provided for that.

## Author

Bishwajit Chakraborty, bishu.math.ynwa@gmail.com
