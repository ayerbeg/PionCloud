# Contributing to PionCloud

PionCloud is developed from several machines by several people (and Claude sessions).
These rules keep its numbers trustworthy.

## 1. Workflow

```bash
git fetch && git status                 # never assume your copy is the latest
git switch -c <topic-branch> origin/main
cmake -S . -B build -DUSE_ROOT=ON -DUSE_CTEQ6=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
```

- One logical change per commit.
- The commit message says **what** and **why**.
- Any user-visible or numerical change gets an entry in `CHANGELOG.md` covering what, why,
  evidence, and whether numbers change.
- Keep `TODO.md` up to date: priority and status, `[open]` or `[done vX.Y]`.
- Pull/rebase before pushing. Never force-push shared branches. Do not commit build products,
  ROOT files or large data.

## 2. The numerical-change policy

Every change falls into one of two classes, and must be shown to belong there.

**(a) Must not change numbers.** Refactoring, comments, I/O, new options with their default
off. Show that all validation runs are **bit-identical** before and after:

```bash
docs/validation/make_runs.sh <old-binary> /tmp/before
docs/validation/make_runs.sh <new-binary> /tmp/after
# compare every numeric column of /tmp/before/runs/*.dat with /tmp/after/runs/*.dat
```

**(b) Changes numbers on purpose** (bug fix, new default, new physics). Document:
- the physical reason, with a reference (paper, equation, table);
- the size of the change on a realistic configuration;
- a unit test that **fails on the old code and passes on the new one**.

Run `docs/validation/run_compare.sh` and check the agreement with the published TDIS
figures (manual, Sec. 7) is unchanged, or explain why it changed.

## 3. Tests

- `tests/test_all.cc` is self-contained (no framework).
- ctest accepts the run **only** if the summary line `Results: N passed, 0 failed` is
  printed. This protects against silent aborts: before v1.4 a Fortran `STOP` ended the
  suite with exit status 0 and ctest reported "Passed".
- Prefer physics tests over "is positive" tests. Good examples:
  - sum rules (`test_PionPDFSumRules`);
  - isospin relations;
  - exactness of quadrature on polynomials;
  - published multiplicities;
  - regression against the Fortran (`test_Regression`).

## 4. Adding physics

- **New splitting function / meson structure function.** Implement `SplittingModel` or
  `PionSFModel` (`include/Models.hh`) and inject or register it. Do not edit the scan loops.
- **New built-in channel.** Follow "HOW TO ADD A CHANNEL" in `include/SplittingFunctions.hh`.
- **Check the normalisation** against a published number. For example:
  ⟨n⟩ = ∫dy ∫dk_T² f(y, k_T) against McKenney et al., PRD 93 054011, Table I
  (`analyses/code_review/multiplicity` on the original laptop shows how).
- **Cut-offs belong to a regulator.** Λ values are not transferable between form factors.
  For example: 1.33 GeV is the s-exponential fit, 0.27 GeV the Pauli–Villars fit.

## 5. Style

- C++17, namespace `PionCloud`, `.hh`/`.cc`.
- Integer-indexed integration loops (`coord = min + i*step`), weights from
  `simpsonWeight()`.
- Units: GeV, and radians internally (degrees only at I/O).
- Comments explain the physics and cite the source (paper and equation, or Fortran file and
  line).
- Block comments at each section say what it does and how to modify it.
