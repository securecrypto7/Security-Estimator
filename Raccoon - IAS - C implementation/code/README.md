# Raccoon-based Interactive Aggregate Signatures (IAS) — reference C implementation

A self-contained, conference-quality reference implementation of the
**Raccoon-based Interactive Aggregate Signature** scheme with **public-key
aggregation** (Fig. 10 of the IAS paper). It implements the full multi-signer
protocol end-to-end, reproduces the security estimator's key/signature sizes
byte-for-byte, and reports per-operation timings for all six parameter sets.

The code has no external dependencies: SHAKE256/Keccak, big-integer ring
arithmetic (up to a ~2^71 NTT prime, held in 128-bit words), the negacyclic
NTT, the Gaussian samplers, the scheme, and the serialization are all included.

---

## Scheme

Ring `R_q = Z_q[X]/(X^n + 1)`, `n = 512`, `q` an NTT-friendly prime
(`q ≡ 1 mod 2n`). `A ←$ R_q^{k×l}` is expanded from a 32-byte seed.

| Round | Party | Output |
|-------|-------|--------|
| `Setup` | — | public matrix `A` (from a seed) |
| `KeyGen` | signer | `sk = s_i`, `pk = t_iᵀ = ⌊2(A s_i + e_i)⌉_{ν_t}` |
| `KeyAgg` | anyone | `a_i = H_agg(PK,i,t_iᵀ) ∈ T`, aggregate key `t̃ = Σ a_i·2^{ν_t}·t_iᵀ` |
| `Sign_pre` | signer | commitments `w_{i,b} = ⌊A r_{i,b}+e'_{i,b}⌉_{ν'_w}`, secret state `r_{i,b}` |
| `Coord` | coordinator | weights `(β_b)=H_β(…)` (β₁=1), aggregate `w = ⌊Σ_j Σ_b 2^{ν'_w}β_b w_{j,b}⌉_{ν_w}` |
| `Sign_resp` | signer | `c = H_c(t̃,(m_i),w)`, `z_i = 2 a_i c s_i + Σ_b β_b r_{i,b}`; partial `σ_i=(c,z_i)` |
| `Coord_combine` | coordinator | per-signer identifiable-abort check, then `z = Σ z_i`, hint `h = w − ⌊A z − c t̃⌉_{ν_w}`; `σ = (c, z, h)` |
| `Verify` | anyone | recompute `w' = ⌊A z − c t̃⌉_{ν_w} + h`, check `c = H_c(t̃,(m_i),w')` and `‖(z, 2^{ν_w} h)‖₂ ≤ B2` |

`H_c` outputs a ternary challenge with `‖c‖∞=1, ‖c‖₁=ω`; `a_i` and `β_b`
are signed monomials `±X^j` (norm-preserving, `β₁` fixed to 1). Secrets are
drawn from rounded Gaussians `D_s` (σ_s = 2^{u_s}) and `D_r` (σ_r = 2^{u_r}).

### Design note — per-signer check uses the *unrounded* key

The interactive per-signer (identifiable-abort) check in `Coord_combine`
computes `y_i = w_i − (A z_i − c a_i t_i)` with the **unrounded** verification
key `t_i = 2(A s_i + e_i)`, exactly as Hermine's partial verification (Fig. 7)
does. If the *rounded* key `2^{ν_t} t_iᵀ` were used here, `y_i` would carry an
extra term `c·a_i·(2^{ν_t} t_iᵀ − t_i)` of magnitude `≈ ω·2^{ν_t−1}·√(nk)`,
which alone exceeds the per-signer bound `B2ind` and would make every honest
run abort. The rounding-induced slack is instead accounted for **once**, at the
aggregate level, in `B2`:

```
B2 = T·B2ind + (3·2^{ν_w−1} + q⊥ + ω·2^{ν_t−1})·√(nk)
```

so `B2ind` deliberately excludes the `t`-rounding term. Each `pk` therefore
stores both `tT` (rounded — the published VK) and `t_full` (unrounded — used
only for the local check). `t_full` is **not** transmitted and does **not**
affect any published size; only the rounded `tᵀ` is aggregated into `t̃`.

---

## Sizes match the security estimator exactly

The verification-key and signature (z-vector) byte counts reproduce
`estimator_FINAL_v2`’s `sizes()` closed forms:

```
VK  = SEED_BYTES + k·n·(⌈log₂ q⌉ − ν_t)/8
SIG = l·n·z_bits/8 ,   z_bits = ⌈log₂(2·B2 + 1)⌉
```

Because `n = 512` is divisible by 8, every polynomial packs into a whole number
of bytes, so the packing is byte-aligned with no padding.

**One subtlety is reproduced deliberately.** The estimator computes
`⌈log₂ q⌉` with Python’s `math.log2`, which first casts `q` to an IEEE-754
double (round-to-nearest). For the two NTT primes that sit only ~2^15 above a
power of two (IAS-192/256 at Q2^60), the double rounds `q` *down* to that power,
so the estimator’s `⌈log₂ q⌉` is one **less** than the true bit length. The
implementation mirrors this in `IAS_LOGQ_SIZE = ⌈log₂((double)q)⌉` for size
accounting, while keeping the exact bit length (`IAS_LOGQ_CEIL`) for uniform
rejection sampling of `A`. This is lossless: the rounded key/hint coefficients
lie in `[0, ⌊q/2^ν⌋)`, which still fits the narrower width (checked by
assertions in the packer).

The estimator’s reported `sig` counts the **z-vector only**. The challenge `c`
and hint `h` are additional on-wire data; the benchmark reports both the
estimator-matching figure and the full on-wire signature.

| Set | z_bits | VK (B) | SIG=z (B) | +c | +h | full sig (B) |
|-----|:------:|:------:|:---------:|:--:|:--:|:------------:|
| IAS-128 Q2^60 | 51 | 10 912 | 19 584 | 24 | 10 560 | 30 168 |
| IAS-128 Q2^30 | 45 |  7 392 | 14 400 | 24 |  7 040 | 21 464 |
| IAS-192 Q2^60 | 56 | 12 320 | 28 672 | 39 | 12 288 | 40 999 |
| IAS-192 Q2^30 | 47 | 11 232 | 21 056 | 39 | 10 752 | 31 847 |
| IAS-256 Q2^60 | 58 | 17 952 | 37 120 | 55 | 14 720 | 51 895 |
| IAS-256 Q2^30 | 49 | 15 584 | 28 224 | 55 | 14 976 | 43 255 |

All rows are asserted equal to the estimator’s `sizes()` at run time
(`VK: … MATCH`, `SIG: … MATCH`), and `z` is verified to pack/unpack bit-exactly.

---

## Benchmarks

`nsigners = 16`, mean of an adaptive timing loop, times in **milliseconds**.
Per-signer operations are reported per single signer (the `×16` column in the
tool output gives the all-signers total). Measured on an Intel Xeon @ 2.10 GHz,
`gcc 13.3.0 -O2 -march=native`, single thread.

| Set | Setup | KeyGen/s | KeyAgg | Sign_pre/s | Coord | Sign_resp/s | Combine | Verify |
|-----|------:|--------:|-------:|-----------:|------:|------------:|--------:|-------:|
| IAS-128 Q2^60 |  5.6 | 1.68 |  81.7 |  23.5 |  337.7 |  8.72 |  36.5 | 2.49 |
| IAS-128 Q2^30 |  4.5 | 1.46 |  76.4 |  20.9 |  335.6 |  7.64 |  31.9 | 2.23 |
| IAS-192 Q2^60 | 13.6 | 3.50 | 136.4 |  77.0 |  985.2 | 23.77 |  92.7 | 5.81 |
| IAS-192 Q2^30 |  8.6 | 2.19 |  99.9 |  46.2 |  612.3 | 14.99 |  49.3 | 3.37 |
| IAS-256 Q2^60 | 21.4 | 4.97 | 188.1 | 133.0 | 1498.3 | 36.02 | 112.3 | 7.38 |
| IAS-256 Q2^30 | 14.7 | 3.03 | 146.0 |  68.0 | 1091.6 | 24.09 |  66.4 | 4.54 |

`KeyGen/s`, `Sign_pre/s`, `Sign_resp/s` are per-signer. `Coord` and `KeyAgg`
scale with `nsigners·rep` and dominate the coordinator cost; `Verify` is a
single de-aggregated NTT check and is by far the cheapest step. These are
straight-line reference numbers (runtime twiddle factors, no AVX, no
precomputation) intended to characterize the protocol’s shape, not to be a
speed record.

Reproduce with `make run` (numbers vary with hardware).

---

## Build & run

Requires a C11 compiler and `libm`. No other dependencies.

```sh
make            # build bench + test for all six parameter sets into build/
make run        # build and run every benchmark (sizes + timings)
make runtest    # build and run every correctness self-test
make clean
```

Build a single set by hand:

```sh
gcc -O2 -march=native -I. -DIAS_PARAM=IAS192_Q60 \
    bench.c ias.c poly.c sample.c serialize.c fips202.c -lm -o bench
./bench
```

`IAS_PARAM` selects the set: `IAS128_Q60`, `IAS128_Q30`, `IAS192_Q60`,
`IAS192_Q30`, `IAS256_Q60`, `IAS256_Q30` (default `IAS128_Q60`).

### What the binaries print

* **`bench_<SET>`** — parameters and bounds, per-operation timings, the
  serialization round-trip check, the size table, and the estimator size
  assertions, ending in `RESULT: PASS`.
* **`test_<SET>`** — a full honest run (`Setup → 16×KeyGen → KeyAgg →
  16×Sign_pre → Coord → 16×Sign_resp → Coord_combine → Verify`), the actual
  aggregate norm vs `B2`, and two negative tests (tampered `z`, wrong message),
  ending in `RESULT: PASS`.

---

## File layout

| File | Contents |
|------|----------|
| `params.h` | all six parameter sets (selected by `-DIAS_PARAM=…`) |
| `fips202.{h,c}` | SHAKE256 / Keccak-f[1600] |
| `poly.{h,c}` | `R_q` arithmetic, 128-bit modular ops, negacyclic NTT, modulus rounding, norms, `q` selection |
| `sample.{h,c}` | seeded XOF RNG, rounded-Gaussian sampler, hashes `H_agg`, `H_β`, `H_c`, matrix-`A` expansion |
| `ias.{h,c}` | the scheme: setup, keygen, keyagg, sign_pre, coord, sign_resp, coord_combine, verify, bounds |
| `serialize.{h,c}` | bit-packing for VK and signature (estimator-matching sizes) |
| `bench.c` | timing + size report + estimator size assertions (main) |
| `test.c` | correctness self-test + negative tests |
| `Makefile` | builds/runs all six sets |

---

## Notes & limitations

This is a **reference** implementation, optimized for clarity and for matching
the paper’s parameters and sizes, not for production. In particular it is not
constant-time, uses `double`-based Gaussian sampling, and computes NTT twiddles
at runtime. The `A`-matrix and all hashes are derived from SHAKE256 with
domain-separated inputs; the RNG for secrets is a seeded Keccak XOF chosen for
reproducible benchmarking and must be replaced by a cryptographic RNG for any
real use.
