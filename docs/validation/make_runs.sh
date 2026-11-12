#!/bin/bash
# Generate PionCloud inputs for every curve of the published reference figures
# (TDIS proposal PR12-15-006, Figs. 5, 6, 37, 46, 47) and run them.
#
# Usage: docs/validation/make_runs.sh [PionCloud binary] [work dir]
#   defaults: <repo>/build/PionCloud   ./validation_runs
# The CTEQ6 .tbl files are copied from the binary's directory into the work dir.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
EXE=$(realpath "${1:-$HERE/../../build/PionCloud}")
WORK=$(mkdir -p "${2:-validation_runs}" && cd "${2:-validation_runs}" && pwd)
IN=$WORK/inputs; RUN=$WORK/runs
mkdir -p "$IN" "$RUN"
cp -n "$(dirname "$EXE")"/cteq6*.tbl "$RUN"/ 2>/dev/null || true

# kbin_method = grid: the ORIGINAL Fortran algorithm, with which the published
# figures were made and the overlays in docs/figures were produced.  The v1.4
# default is exact (see CHANGELOG.md); tests/validation.cc checks both.
COMMON="kbin_method = grid
form_factor = exp_s
lambda = 1.56
lambda_delta = 1.39
beam_energy = 11.0
output_mode = all
verbose = false"

# mk name key=val ...   -> writes inputs/name.in
mk() { local n=$1; shift; { echo "$COMMON"; echo "output_name = $n"; for kv in "$@"; do echo "${kv/=/ = }"; done; } > "$IN/$n.in"; }

# ---- Fig. 47 (PR12-15-006): theta_p' dependence, neutral exchange, theta_e = 35 deg (3Var_theta.f)
for ch in pion rho; do
  mk fig47L_$ch scan_mode=theta contribution=$ch dissociation=neutral electron_angle=35 \
     theta_min=0 theta_max=90 nth=90 n_x_integration=200 nk=200
  mk fig47R_$ch scan_mode=theta contribution=$ch dissociation=neutral electron_angle=35 \
     theta_min=0 theta_max=90 nth=90 n_x_integration=200 nk=200 \
     x_min_th=0.05 x_max_th=0.6 k_min_th=0.060 k_max_th=0.250
done

# ---- Fig. 46: x dependence, Delta|k| = [0,500] MeV, all angles, neutral (3Var_x.f)
for ch in pion rho; do
  mk fig46_$ch scan_mode=x contribution=$ch dissociation=neutral electron_angle=35 \
     x_min=0.005 x_max=0.4 nx=80 theta_min=0 theta_max=180 n_theta_integration=200 \
     k_min=0.0001 k_max=0.5 nk=200
done

# ---- Fig. 5 (= Hobbs 1411.2216 Fig. 2 right): p -> pi0 p, four |k| bins, F2p
#      two angular windows: header [30,180] deg and the 30-70 deg window of the repo inputs
for kb in "060 0.060 0.100" "100 0.100 0.200" "200 0.200 0.300" "300 0.300 0.400"; do
  set -- $kb
  mk fig5_k$1_th30-180 scan_mode=kbin contribution=pion dissociation=neutral nucleon=proton \
     electron_angle=35 x_min=0.05 x_max=0.30 nx=50 k_min=$2 k_max=$3 \
     cosph_min=0.0 cosph_max=0.866 ny=100 nkt=5000
  mk fig5_k$1_th30-70 scan_mode=kbin contribution=pion dissociation=neutral nucleon=proton \
     electron_angle=35 x_min=0.05 x_max=0.30 nx=50 k_min=$2 k_max=$3 \
     cosph_min=0.342 cosph_max=0.866 ny=100 nkt=5000
done
mk fig5_inclusive scan_mode=kbin contribution=pion dissociation=neutral nucleon=proton \
   electron_angle=35 x_min=0.05 x_max=0.30 nx=50 inclusive_mode=true ny=100 nkt=5000

# (cos(phi_k) > 0 always in this kinematics, so [30,180] deg == cosph in [0, 0.866])
# ---- Fig. 6: n -> pi- p, rho- p, pi Delta; Delta theta = [30,180] deg; Delta|k| = [60,250] MeV
# ---- Fig. 37 (= Hobbs Fig. 2 left): same channels, Delta|k| = [250,400] MeV, Delta theta = [30,70] deg
for ch in pion rho pi_delta; do
  mk fig6_$ch scan_mode=kbin contribution=$ch dissociation=charge_exchange nucleon=neutron \
     electron_angle=35 x_min=0.05 x_max=0.20 nx=50 k_min=0.060 k_max=0.250 \
     cosph_min=0.0 cosph_max=0.866 ny=100 nkt=5000
  mk fig37_$ch scan_mode=kbin contribution=$ch dissociation=charge_exchange nucleon=neutron \
     electron_angle=35 x_min=0.05 x_max=0.30 nx=50 k_min=0.250 k_max=0.400 \
     cosph_min=0.342 cosph_max=0.866 ny=100 nkt=5000
done
for fig in fig6 fig37; do
  mk ${fig}_inclusive scan_mode=kbin contribution=pion dissociation=charge_exchange nucleon=neutron \
     electron_angle=35 x_min=0.05 x_max=0.30 nx=50 inclusive_mode=true ny=100 nkt=5000
done

cd "$RUN"
for f in "$IN"/*.in; do
  n=$(basename "$f" .in)
  "$EXE" --input "$f" > "log_$n.txt" 2>&1 || { echo "FAILED: $n"; tail -3 "log_$n.txt"; }
done
echo "done: $(ls "$RUN"/*_xscan.dat "$RUN"/*_thetascan.dat 2>/dev/null | wc -l) data files"
