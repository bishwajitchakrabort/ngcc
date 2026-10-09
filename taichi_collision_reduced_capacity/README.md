# TaiChi IFS structural collision attack (real P1920, reduced capacity)


## Disclaimer

The collision attack described in the accompanying paper was produced by
the author. It was subsequently verified with the assistance of an AI system,
and all verification code (the reimplementation of both WChain versions and
their inverse permutations, the comparison with the reference implementation,
the reduced-width attack, and the ideal-permutation experiment) was written by
that AI system. The paper was itself drafted in its entirety by that AI system.

## Attack

Structural 2^(3c/4) collision attack on the TaiChi IFS mode, run against the **real**
TaiChi permutation `P1920` (P1 = round-constant bank 1, P2 = bank 2), with the capacity
reduced to `c` bits and the rate set to `1920 - c` bits as in any sponge. The permutation is
the unmodified NGCC reference; only the mode-level rate/capacity split is scaled down so the
attack is runnable.

This is not a generic birthday search. It is the meet-in-the-middle + birthday attack on the
cumulative-XOR relation R1 = x^y, final P2 capacity input = x^y^z = D: fix the rate part A*
with inverse P1 queries, collide the c-bit D. The real TaiChi-512 has c = 512, so the attack
is 2^384 there and does not threaten it; these runs verify the scaling law and produce
concrete collisions at small c.

## Files

- `perm.hpp`            real P1920, forward + inverse, scalar and batched (BS=64 lanes)
- `taichi_real_attack.cpp`  mode, attack, engines, selftest, bench, scale, batch
- `verify.py`          independent Python reimplementation; checks a collision from scratch
- `ref.c`              the NGCC reference (verbatim) + a driver, used as ground truth
- `build.sh`           build helper (Apple clang / Homebrew gcc / Linux)

## Build

```
sh build.sh                       # auto flags
# or explicitly on macOS:
clang++ -O3 -std=c++17 -pthread taichi_real_attack.cpp -o taichi_real
# with Homebrew gcc (wider SIMD, usually faster):
g++-14 -O3 -mcpu=native -std=c++17 -pthread taichi_real_attack.cpp -o taichi_real
```

`c` must be a multiple of 8 in [8, 64]. The code uses all cores by default (`--threads N`
to override).

## First: convince yourself it is correct

```
./taichi_real selftest
```

Checks, in order: P1920 forward matches the reference vectors; the mode reproduces the four
**real TaiChi-512** reference digests (so the permutation and the mode wiring are exact);
inverse round-trips; the reduced mode's first-block state is (s, y, x^y); 2R-2 bit messages
pad back to the two chosen blocks; the partitioned engine equals brute force and agrees with
the interleaved engine on the first-collision q; a full c=16 collision verified through the
reference hash.

Confirm the exponent on the real permutation (a few minutes):

```
./taichi_real scale --cs 8,16,24
```

Prints median/mean log2 q per c and a weighted fit of the exponent (theory 0.75).

## Run c = 40

```
./taichi_real batch --c 40 --mem-gb 11 --out c40
```

Set `--mem-gb` to about 60-70% of your RAM. The planner prints the plan, the balanced-optimum
vs this-plan query counts, and a **measured ETA** from a short on-machine calibration before it
starts. It makes one pass when `q_f` fits in the table, otherwise splits into passes. On the
first D-collision it stops, prints both messages (and writes `c40_M.hex`, `c40_Mprime.hex`),
recomputes both hashes through the reference mode, and asserts the full 1920-c bit final state
(hence the digest) is identical.

Independent check:

```
python3 verify.py 40 c40_M.hex c40_Mprime.hex
```

Reimplements P1920 and the mode in Python (no shared code) and confirms the digests collide.

Useful switches: `--lambda L` expected D-collisions (default 4, P[success] 98%; raise for more
certainty, lower for speed), `--seed S`, `--passes P`, `--no-stop` (enumerate all survivors),
`--threads N`, `--q Q` / `--qf Q --qi Q` (set the query counts by hand).

## Parallelism

The batch engine runs the forward table-build and the inverse stream across all threads. The
D-map is split into 256 stripes (equal D always shares a stripe, so detection stays exact) and
each thread keeps a private survivor buffer merged at the end, so the only shared-state traffic
is the rare survivor insert (about 1 in 10^4 queries). `bench` reports throughput on 1 and all
threads.

## Feasibility (why c = 40 is fine and c = 64 is not)

Balanced cost at target lambda = 4 is q = 2^(0.75c + 1.75) forward roots and the same number of
inverse queries, i.e. ~3 * 2^(0.75c + 1.75) permutation calls, and the join table wants 8 * q
bytes. Memory is what bites first; beyond the one-pass limit the planner trades it for time by
making q_i larger than q_f (more inverse queries), or you add `--passes`.

| c  | balanced q/side | total P1920 calls | balanced table | on an M3 Air (8 cores) |
|----|-----------------|-------------------|----------------|------------------------|
| 32 | 2^25.8          | 2^27.3  (~1.7e8)  | 0.45 GB        | seconds                |
| 40 | 2^31.8          | 2^33.3  (~1.1e10) | 29 GB          | ~10-25 min (1 pass, memory-limited q_i) |
| 44 | 2^34.8          | 2^36.3  (~8.6e10) | 230 GB         | a few hours (overnight), memory-limited |
| 48 | 2^37.8          | 2^39.3  (~6.9e11) | 1.8 TB         | ~weeks, not practical  |
| 64 | 2^49.8          | 2^51.3  (~2.8e15) | 7.6 PB         | no: see below          |

**c = 64 is not possible on a single machine.** It needs ~2^49-2^50 queries per side (about
10^15 permutation evaluations ~ order a year on one laptop even at full tilt) and, balanced,
petabytes of join table. Keeping q_i <= 2^c (the only way to bound memory to a laptop) forces
q_f >= ~2^35.5 entries = ~384 GB of table and q_i up to 2^64 inverse queries, which is
millennia. The attack's whole point is this 2^(3c/4) growth: each +4 in c multiplies work by 8
and balanced memory by 16. At the real c = 512 it is 2^384.

Rule of thumb on an M3 Air: c <= 40 is comfortable, c = 44 is an overnight run, c >= 48 is out
of reach, and run `./taichi_real bench --c 40` first — the batch planner's ETA line is computed
from exactly that measurement on your hardware.
