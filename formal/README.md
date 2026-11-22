# Lean check of the high-zero integral

The [Lean file](TailIntegration.lean) formalizes the repeated integration by parts in the [high-zero estimate](../tapered-bound.md#controlling-zeros-above-the-verified-height). Given a continuous chain of derivatives and a function that vanishes outside the integration interval, it derives the vanishing endpoint terms and proves the bound for a single zero. With 390 integrations by parts, the denominator has the 391st power of the absolute imaginary ordinate.

The check takes the derivative chain and its final bound as inputs. The construction and derivative estimates for the specific taper, the sum over all high zeros, and the crossover conclusion are documented in the research note.

## Reproduce

Use Lean 4.28.0 and Mathlib v4.28.0. From the root of this repository, run the following commands.

```sh
mkdir -p data/lean-check
git clone --depth 1 --branch v4.28.0 https://github.com/leanprover-community/mathlib4.git data/lean-check/mathlib
cd data/lean-check/mathlib
lake exe cache get Mathlib.Analysis.SpecialFunctions.ExpDeriv Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts Mathlib.Tactic.Ring
lake env lean -DwarningAsError=true ../../../formal/TailIntegration.lean
```

The final command prints the axioms used by each theorem.
