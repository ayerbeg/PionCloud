# PionCloud — Pion Cloud Model Structure Function Calculator

C++ port of `3Var_x.f` and `3Var_theta.f` (T. Hobbs / W. Melnitchouk).
Computes the leading pion-cloud contribution to the nucleon structure
function F₂ in two complementary scan modes:

| Mode | Mirrors | Outer variable | Integrated over |
|------|---------|----------------|-----------------|
| `x`     | `3Var_x.f`     | Bjorken-x  | θ_h and \|k\| |
| `theta` | `3Var_theta.f` | θ_h (deg)  | x and \|k\|   |
| `all`   | both           | both       | —              |

---

## Project layout

```
PionCloud/
├── CMakeLists.txt
├── CTEQ6/
│   ├── Cteq6Pdf-2007.f   ← Fortran PDF library source (place here)
│   ├── cteq6m.tbl        ← grid table for CTEQ6M (iset=1)
│   ├── cteq6d.tbl        ← grid table for CTEQ6D (iset=2)
│   ├── cteq6l.tbl        ← grid table for CTEQ6L (iset=3)
│   ├── cteq6l1.tbl       ← grid table for CTEQ6L1 (iset=4)
│   └── README            ← setup instructions
├── include/
│   ├── PhysicsParams.hh      ← masses, enums, RunParams struct
│   ├── SplittingFunctions.hh ← fypiN, f_rhoN, f_RhoDel, fypiD (inline)
│   ├── GRV.hh                ← GRV92 LO pion PDFs (inline)
│   ├── Kinematics.hh         ← Q2, (k,θ)→(y,kT), phase space, Simpson weights
│   ├── Results.hh            ← XScanPoint, ThetaScanPoint, ScanResults
│   ├── Calculator.hh         ← integration engine (interface + scan constants)
│   ├── CTEQ6.hh              ← extern C bridge to Fortran library
│   ├── InputReader.hh        ← INI file parser (no physics dependency)
│   └── OutputWriter.hh       ← ASCII + ROOT writer
└── src/
    ├── Calculator.cc         ← integer-indexed Simpson integration loops
    ├── InputReader.cc        ← parser, validator, template writer
    ├── OutputWriter.cc       ← file writers
    └── main.cc               ← CLI entry point
```

---

## Building

### Prerequisites
- CMake ≥ 3.14  **or**  g++ ≥ 7 directly
- C++17 or later
- gfortran (for CTEQ6; optional but recommended)
- ROOT ≥ 6 (optional; needed only for `.root` output)

### Full build (ROOT + CTEQ6)
```bash
source /path/to/root/bin/thisroot.sh
cmake -B build -DUSE_ROOT=ON -DUSE_CTEQ6=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
# binary: build/PionCloud
# tests:  build/PionCloudTests
```

### Without ROOT
```bash
cmake -B build -DUSE_ROOT=OFF -DUSE_CTEQ6=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

### Without CTEQ6 (F2n and ratio will be zero)
```bash
cmake -B build -DUSE_ROOT=ON -DUSE_CTEQ6=OFF
cmake --build build
```

### Direct g++ (no CMake, no CTEQ6, no ROOT)
```bash
g++ -std=c++17 -O3 -Iinclude \
    src/Calculator.cc src/InputReader.cc src/OutputWriter.cc src/main.cc \
    -o PionCloud
```

With ROOT and CTEQ6 via g++ directly:
```bash
gfortran -O2 -c CTEQ6/Cteq6Pdf-2007.f -o Cteq6Pdf.o
g++ -std=c++17 -O3 -Iinclude -DHAVE_ROOT -DHAVE_CTEQ6 \
    src/Calculator.cc src/InputReader.cc src/OutputWriter.cc src/main.cc \
    Cteq6Pdf.o \
    $(root-config --cflags --libs) -lgfortran \
    -o PionCloud
```

### Running tests
```bash
ctest --test-dir build -V
# or directly:
./build/PionCloudTests
```

---

## Running

### Step 1 — generate a template input file
```bash
./PionCloud --write-template myrun.in
```

### Step 2 — edit the input file
Every key is optional; unset keys keep their defaults.
Example for a rho-N x-scan with BoNuS kinematics:

```ini
[Run Control]
scan_mode    = x

[Physics - Form Factor]
form_factor  = exp_s        # typ=2

[Physics - Meson-Baryon Channel]
contribution = rho

[Physics - Cutoff Parameters]
lambda       = 1.56         # HSS central value

[x-Scan Settings]
x_min = 0.02
x_max = 0.30
nx    = 100

[|k| Integration]
k_min = 0.060
k_max = 0.160
nk    = 100
```

### Step 3 — run
```bash
./PionCloud --input myrun.in
```

### CLI flags override the input file
```bash
# Sweep over BoNuS |k| bins
for k1 in 0.060 0.080 0.100 0.130 0.160; do
    k2=$(echo "$k1 + 0.020" | bc)
    ./PionCloud --input myrun.in --kmin $k1 --kmax $k2 --out "Run_k${k1}"
done
```

---

## All options

### Input file / template
```
--input  <file>              Read parameters from INI-style file
--write-template <file>      Write a documented template and exit
```

### Run control
```
--mode  x|theta|all          Scan mode (default: all)
```

### Output
```
--output  ascii|root|all|none    Format (default: all)
--out     <basename>             File stem (default: PionCloud)
```

### Physics
```
--flag   0|1|2     0=pi-N, 1=rho-N, 2=pi-Delta
--dis    0|1       0=charge exchange, 1=neutral
--typ    0..6      Form factor (see table below)
--L    <GeV>       Λ for pi/rho-N vertex
--Ld   <GeV>       Λ for Delta vertex
--E    <GeV>       Beam energy
--theta_e  <deg>   Electron scattering angle
--pH   <GeV>       Produced hadron momentum
--alpha1   <deg>   Hadron angle lower bound
--alpha2   <deg>   Hadron angle upper bound
```

### Form factor types
| Code | Name | Notes |
|------|------|-------|
| 0 | `monopole` | |
| 1 | `dipole` | |
| 2 | `exp_s` | s-dep. exponential; HSS convention |
| 3 | `cov_dipole` | Covariant dipole |
| 4 | `dipole_sch` | Dipole, s-channel Λ exchange |
| 5 | `exp_t` | t-dep. exponential; PRD93,054011(2016) |
| 6 | `pauli_villar` | Pauli-Villar; PRD93,054011(2016) **[default]** |

### x-scan
```
--xmin <val>        Minimum Bjorken-x  (default: 0.02)
--xmax <val>        Maximum Bjorken-x  (default: 0.30)
--nx   <n>          Output x points    (default: 100)
--nth-int <n>       Internal θ_h steps (default: 100)
```

### theta-scan
```
--thmin <deg>       θ_h lower bound     (default: 0.0)
--thmax <deg>       θ_h upper bound     (default: 100.0)
--nth   <n>         Output θ points     (default: 100)
--nx-int  <n>       Internal x steps    (default: 100)
```
> The theta-scan integrates x over [~0, 0.6] and k over [~0, 0.5 GeV]
> (full phase space, as in `3Var_theta.f`). `nk` and `nx-int` control resolution.

### |k| integration
```
--kmin <GeV>        Lower |k| bound (x-scan tagging window, default: 0.05)
--kmax <GeV>        Upper |k| bound (x-scan tagging window, default: 0.10)
--nk   <n>          Steps, used in both scans (default: 100)
```

---

## Output files

With `--out MyRun`:

| File | Content |
|------|---------|
| `MyRun_params.txt`    | Full parameter record + first moments |
| `MyRun_xscan.dat`     | `x  Q2  F2piK  RunningMoment  [F2n  ratio]` |
| `MyRun_thetascan.dat` | `theta_h_deg  F2piK  RunningMoment` |
| `MyRun.root`          | TGraphs, TH1s, TCanvases, TNamed params |

Columns in `[]` appear only when built with CTEQ6.

### ROOT objects

Object names encode the **channel and scan mode** so multiple ROOT files
can be merged without name collisions. The prefix is `<channel>_<mode>`:
- channel: `piN` | `piN_neutral` | `rhoN` | `piDelta`
- mode: `xscan` | `theta` | `kbin`

**Always present when xScan non-empty:**

| Name | Type | Description |
|------|------|-------------|
| `g_<pfx>` | `TGraph` | F₂^channel(x), e.g. `g_piN_kbin` |
| `g_moment_<pfx>` | `TGraph` | Running ∫F₂ dx |
| `h_<pfx>` | `TH1D` | Same as histogram |
| `c_<pfx>` | `TCanvas` | Log-scale plot |

**Present when built with CTEQ6:**

| Name | Type | Description |
|------|------|-------------|
| `g_F2neutron` | `TGraph` | F₂ⁿ(x) from CTEQ6M — **always this name** |
| `g_ratio_<pfx>` | `TGraph` | F₂^channel / F₂ⁿ, e.g. `g_ratio_piN_kbin` |
| `h_ratio_<pfx>` | `TH1D` | Ratio as histogram |
| `c_ratio_<pfx>` | `TCanvas` | Two-pad: curves (top) + ratio (bottom) |

**Present when thetaScan non-empty:**

| Name | Type | Description |
|------|------|-------------|
| `g_thetascan` | `TGraph` | F₂^π(θ_h) |
| `g_moment_theta` | `TGraph` | Running ∫F₂ dθ |
| `h_thetascan` | `TH1D` | Same as histogram |
| `c_thetascan` | `TCanvas` | Angular distribution plot |
| `RunParams` | `TNamed` | Full parameter summary string |

**Which name corresponds to which curve in the f2piN plot:**

| Curve | ROOT object name |
|---|---|
| Orange F₂ⁿ(x) | `g_F2neutron` |
| Purple F₂^(πN)(x) inclusive | `g_piN_kbin` (no angle/k cuts) |
| Black πN with acceptance | `g_piN_kbin` (30°–70°, |k|=[60,250] MeV) |
| Red ρN | `g_rhoN_kbin` |
| Green πΔ | `g_piDelta_kbin` |
| Fig.47 θ distribution | `g_thetascan` |

---

## CTEQ6 — neutron F₂ and ratio F₂^π / F₂ⁿ

### Setup

Place the following files in `CTEQ6/` (all four are part of your distribution):

```
CTEQ6/
  Cteq6Pdf-2007.f   ← Fortran source (v6.52, April 2007)
  cteq6m.tbl        ← CTEQ6M  NLO MSbar  iset=1  ← used by original Fortran
  cteq6d.tbl        ← CTEQ6D  NLO DIS    iset=2
  cteq6l.tbl        ← CTEQ6L  LO         iset=3
  cteq6l1.tbl       ← CTEQ6L1 LO (αs fitted) iset=4
```

See `CTEQ6/README` for the full iset table and compiler notes.

### Physics

F₂ⁿ is computed from the **proton** CTEQ6M PDFs via **charge symmetry**
(u_neutron = d_proton, d_neutron = u_proton), matching `f2_pi_sub.f90` line 269:

```
F2neu = x * ( (4/9)*(d_pro + dbar_pro) + (1/9)*(u_pro + ubar_pro) )
```

The default PDF set is CTEQ6M (iset=1, NLO MSbar), matching the original
Fortran `CALL SETCTQ6(1)  ! CTEQ 'MS-bar' SCHEME`. This can be changed
programmatically via `CTEQ6::init(CTEQ6::CTEQ6L)` etc.

### Runtime note

`SetCtq6` opens the table file by name in the **current working directory**.
CMake copies the `.tbl` files to the build directory at configure time, so
running from `build/` works automatically. If you run the executable from
a different directory, symlink or copy the `.tbl` files there.

---

## Using as a library in another program

```cmake
# In your CMakeLists.txt:
add_subdirectory(PionCloud)
target_link_libraries(MyProgram PRIVATE PionCloudLib)
```

```cpp
#include "PhysicsParams.hh"
#include "InputReader.hh"
#include "Calculator.hh"
#include "Results.hh"

using namespace PionCloud;

// Load from input file
InputReader reader;
RunParams p = reader.read("myrun.in");

// Or configure programmatically
RunParams p;
p.flag       = FLAG_RHO;
p.typ        = FF_EXP_S;
p.L          = 1.56;
p.kmin       = 0.060;
p.kmax       = 0.100;
p.outputMode = OUT_NONE;   // no files written

// Run and access results directly
Calculator calc(p);
ScanResults res;
calc.runXScan(res);

for (auto &pt : res.xScan)
    printf("x=%.4f  F2piK=%.6e  F2n=%.6e  ratio=%.6e\n",
           pt.x, pt.F2piK, pt.F2n, pt.ratio);
```

---

## Physics notes

All splitting functions are translated directly from the Fortran originals:

| Function | Description |
|----------|-------------|
| `fypiN`    | π-N; TOPT; Melnitchouk (1999) / Hobbs (2013) |
| `f_rhoN`   | ρ-N; TOPT with P−p derivative coupling |
| `f_RhoDel` | ρ-Δ; TOPT |
| `fypiD`    | π-Δ; SU(4) analogue; Holzenkamp et al. |
| `GRV92`    | LO pion PDFs; Gluck, Reya, Vogt, Z.Phys.C53(1992)651 |

The coordinate change (k, θ_h) → (y, kT) and the phase-space Jacobian
follow the Fortran exactly. All integrations use the composite Simpson rule
with integer-indexed loops to avoid floating-point loop-variable accumulation.


---
## DISCLOSURE

This project used Anthropic Claude to convert the fortran code into a C++ code. For the moment, it has been tested reproducing several plots already produced for several people. The code is still under development to get more documentation and flexibility to modify as needed. The ONLY AND UNIQUE authors of the algorithm are the original coders of the fortran programs. 