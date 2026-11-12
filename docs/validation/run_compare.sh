#!/bin/bash
# Convenience wrapper: build inputs, run PionCloud and produce the overlay plots.
# Usage: docs/validation/run_compare.sh [PionCloud binary] [work dir]
# Needs PyROOT and numpy (source thisroot.sh first if ROOT is not in the environment).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${2:-validation_runs}
"$HERE/make_runs.sh" "${1:-$HERE/../../build/PionCloud}" "$WORK"
python3 "$HERE/compare.py" "$WORK" | grep -v "^Info in"
