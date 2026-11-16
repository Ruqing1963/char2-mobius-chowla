# char2-mobius-chowla

Paper, verification code and figures for

> Ruqing Chen, *Short-interval Möbius sums, Chowla correlations, and polynomial Möbius sums in
> characteristic 2* (2026).

Archived at DOI [10.5281/zenodo.24267751](https://doi.org/10.5281/zenodo.24267751).
This is the sequel to [short-intervals-char2](https://github.com/Ruqing1963/short-intervals-char2)
(DOI [10.5281/zenodo.24222056](https://doi.org/10.5281/zenodo.24222056)), whose parity law it uses.

## Setting

Let `q = 2^k` and let `I_2(f) = f + {polynomials of degree <= 2}` be a short interval around a monic `f`
of degree `n = 2d >= 4` over `F_q` with `deg f' <= 1`. Every separable member is `F = H^2 + a t`, and by
the parity law of the first paper

    mu(H^2 + a t) = (-1)^(k*eps_d) * chi(tau(H)/a),    tau(H) = sum over roots beta of H of beta*H'(beta)^2,

where `chi(x) = (-1)^Tr(x)` and `eps_d = 1` iff `d = 2, 3 (mod 4)`. Write `tau_f(B0, B1)` for `tau` on
the plane `H = H_f + B1*t + B0`.

## What the paper proves

* **Structure of `tau_f`.** For `d` odd it is bi-affine: `tau_f = (B0 + c0)(B1 + c1) + Delta`. For `d`
  even it is `F_2`-affine: `tau_f = T01*B0 + T1*B1 + h_(d-1)*B1^2 + T00` (from the first paper).
* **Exact Möbius sums.** The sum of `mu` over `I_2(f)` is `(-1)^(k*eps_d) * q * (N(tau_f) - q)`, where `N`
  counts the zeros of `tau_f`. It lies in `{0, ±q, ±q(q-1), ±q^2, ±q^2(q-1)}` and is at most `q^2` in
  absolute value when `tau_f != 0`. Exact formulas are also given on `f + F_q` and
  `f + {degree <= 1}`.
* **Chowla correlations on `I_2(f)`.** An exact dichotomy between `O(q^(5/2))` cancellation and full
  correlation `±q^2(q - |A|)`.
  * There are intervals with geometric Galois group `S_n` on which `mu(F + c) = mu(F)` for every
    constant `c`, so 2-point sums with constant shifts do not cancel.
  * For every `n = 0 (mod 4)`, `n >= 8`, 4-point sums over parallelogram shifts do not cancel.
  * For `d` even the sums are evaluated exactly, without the Weil bound: cancellation is
    all-or-nothing (Remark 3.2).
* **Polynomial Möbius sums.** `|sum over monic H of degree d of mu(H^2 + a t)| <= (d-1) q^(d-1)` for
  `d >= 3` (`<= q^(d-1)` for `d` odd), with closed formulas for `d <= 6`; for example the sum is
  `(-1)^k q^4 + q^3(q-1)` for `d = 6`. This is a characteristic-2 counterpart of a problem of Sawin and
  Shusterman.
* **Rigid intervals.** For every `n = 8m >= 16`, a `2m`-parameter family on which `tau_f` is constant;
  the rigid intervals form a hypersurface in it. For `n = 16` this is the sporadic rigid component
  of the first paper.

## Contents

| path | content |
|---|---|
| `mobius_chowla_char2.tex`, `.pdf` | the paper (11 pages) |
| `figures/` | the three figures (PDF) |
| `code/verify_mobius_chowla.py` | exhaustive checks of the main theorems by direct factorisation (pure Python) |
| `code/extra_checks.py` | supplementary checks, sections 1 to 11 (pure Python; takes several minutes) |
| `code/make_figures.py` | regenerates the figures (numpy, matplotlib) |

## Reproducing

```
python code/verify_mobius_chowla.py   # no dependencies
python code/extra_checks.py           # no dependencies
python code/make_figures.py           # optional: writes figures/*.pdf
pdflatex mobius_chowla_char2.tex      # twice
```

`extra_checks.py` covers, by direct factorisation over `F_2`, `F_4`, `F_8` and `F_16`:

| section | checks |
|---|---|
| 1 | the parity law |
| 2 | the case formulas of Theorem 1.2 (a)-(c) |
| 3, 4a, 4b | Theorem 1.3 (a)-(c): random shift sets, the `n = 8m+4` family, constant shifts |
| 4c | Remark 2.2: the rigid branch for `n = 16` and its count `(q-1)q^2` |
| 5 | the cubic sums `sum_b chi(b^3/a)` |
| 6, 9 | Theorem 1.4: the formulas for `d` odd and even, and the closed forms for `d = 5, 6` |
| 7 | Theorem 1.2(d): the very short intervals |
| 8 | Theorem 1.3(d): parallelogram shifts |
| 10 | Proposition 2.1: the rigid hypersurface for `n = 16` to `40` |
| 11 | Remark 3.2: exact evaluation and bounds of Chowla sums for `d` even |

## Figures

1. `fig_planes.pdf`: the zero locus of `tau_f` and the values of `mu` on a plane over `F_16`, in the
   four typical cases.
2. `fig_chowla.pdf`: normalised Chowla correlations for `q = 4` to `128`, with and without cancellation.
3. `fig_polysums.pdf`: the closed forms for the polynomial Möbius sums, `d = 2` to `6`.

## Citation

```
@misc{Chen2026MobiusChowla,
  author = {Ruqing Chen},
  title  = {Short-interval M\"obius sums, Chowla correlations, and polynomial M\"obius sums in characteristic 2},
  year   = {2026},
  doi    = {10.5281/zenodo.24267751}
}
```
