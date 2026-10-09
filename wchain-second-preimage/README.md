# WChain second-preimage attack (toy verification)

## Disclaimer

The second-preimage attack described in the accompanying paper was produced by
the author. It was subsequently verified with the assistance of an AI system,
and all verification code (the reimplementation of both WChain versions and
their inverse permutations, the comparison with the reference implementation,
the reduced-width attack, and the ideal-permutation experiment) was written by
that AI system. The paper was itself drafted in its entirety by that AI system.

## Attack
A fixed-prefix second-preimage attack on **WChain** (NGCC candidate hash;
Version I `b=576`, `h=512`; Version II `b=1152`, `h=1024`). A four-block
palindrome prefix `X || Y || Y || X` with `P_X(0) = P_Y^{-1}(0)` makes
`Hash(X||Y||Y||X || M) = Hash(M)` for every message `M`. One pair `(X,Y)` is
found in about `2^(b/2)` calls and serves every target, giving second preimages
at about `2^288` (Version I) and `2^576` (Version II), below the `2^512` and
`2^1024` second-preimage requirements. The analysis is in the accompanying
paper.

This folder holds a toy-scale verification. It runs the exact WChain mode (the
two-step feedback, the checksum block `Sigma`, and `10*` padding) at reduced
state width with a generic permutation, so the prefix attack can be checked on a
laptop in under a second.

## Contents

- `README.md` - this file
- `wchain_2ndpreimage_toy.py` - the verification script (no dependencies)

## Requirements

Python 3.8 or later. Standard library only.

## Run

```
python3 wchain_2ndpreimage_toy.py
```

Edit `B_BITS` at the top of the script to change the toy state width.

## What it checks

Sample output at `b = 20`:

```
IV=0  XYYX     : search calls = 2^11.1 (b/2 = 10), Hash(XYYX||M)==Hash(M) on all trials: True
IV!=0 XYYX only: Hash equal on all trials: False (expected False)
IV!=0 XYYX ZWWZ: search calls = 4339 (vs 2162 for one match), Hash equal on all trials: True
```

- **IV = 0**: a single palindrome `XYYX` holds. One pair is found in about
  `2^(b/2+1)` calls, and `Hash(XYYX||M) = Hash(M)` on every random message.
- **IV != 0, XYYX only**: a single palindrome only swaps the register pair, so
  it does not hold (printed as the expected `False`).
- **IV != 0, XYYX then ZWWZ**: a second palindrome closes the loop, at twice the
  cost.



