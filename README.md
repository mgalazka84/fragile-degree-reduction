# Exact verification for degree reduction in fragile power domination

This repository contains the verification code accompanying the manuscript **Degree reduction and heterogeneous
failures in power domination**, by Marek Gałązka and Hanna Wdowicka.

## Run

Python 3.10 or later is sufficient. No external packages or datasets are needed.

```sh
python3 verify.py --output verification.json
```

On Windows, or when the Python 3 interpreter is named `python`, use:

```sh
python verify.py --output verification.json
```

The script raises `AssertionError` on a mismatch and prints a summary on success.
It uses exact integers and `fractions.Fraction`; no random sampling is performed.

## Files

- `fragile.py`: power-domination closure, Möbius coefficients, layer transforms,
  the closed alpha coefficients, elementary tree components and the compressed
  forest recipe and its factorial size bound. `forest_recipe(d)` uses gcd-reduced
  cancellation multipliers; `forest_recipe(d, reduce_multipliers=False)`
  reproduces the original recipe.
- `verify.py`: finite checks, including an independently written set-based
  propagation routine, direct probability enumeration and coefficient checks.
- `verification.json`: the complete output from the delivered version.

## Scope

The enumeration includes every labelled simple graph on 1 through 5 vertices and
every sensor subset, including the empty placement and disconnected graphs. It
does not enumerate all graphs on larger orders. The 18-vertex example is fully
enumerated over its 256 sensor survival patterns. Larger forest constructions
are verified through their compressed component multiplicities, not through
enumeration of their exponentially many survival patterns.

The checks also compare the closed coefficients with rescaled Lagrange weights
for every admissible parameter tuple with `s <= 14`. Both forest recipes are
checked against the proved upper bound `(5/2) d! (d+1)!` for `2 <= d <= 30`.
The factorial-weighted coefficient estimates used in that proof are checked
for both tree families through `j = 100`. These finite ranges supplement the
mathematical proof; they do not establish the bound for all degrees by testing.

The finite checks are supporting evidence, not a formal proof of the theorems.
The code is intended for reproduction and inspection, not for large-scale sensor
placement optimization. The manuscript and its LaTeX source are distributed separately from this code repository.

## Version used for the article

The manuscript identifies the exact Git commit used for its reported checks.
The saved `verification.json` contains the deterministic output of that code.

## License

The code and accompanying repository documentation are available under the
MIT License; see `LICENSE`.
