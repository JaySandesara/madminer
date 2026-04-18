#!/bin/bash
set -e

# First arg (if given) is HTCondor's $(Process) — converted to 1-indexed BATCH_ID.
# Default to 1 when run locally without args.
if [ -n "${1:-}" ]; then
    export BATCH_ID=$(($1 + 1))
else
    export BATCH_ID=${BATCH_ID:-1}
fi
echo "BATCH_ID=$BATCH_ID"

export SETUPTOOLS_SCM_PRETEND_VERSION=0.0.0
python3 -m pip install --no-deps --user -e .

echo "PATH=$PATH"
command -v python3
python3 -c "import sys; print(sys.executable)"

# # madminer is already installed in the container image (same version 0.6.3.dev668)
# python3 -c "import madminer; print('madminer location:', madminer.__file__)"
# echo ""

if [ -d tutorial_particle_physics ]; then
    cd tutorial_particle_physics
elif [ -d examples/tutorial_particle_physics ]; then
    cd examples/tutorial_particle_physics
else
    echo "ERROR: cannot find tutorial_particle_physics directory" >&2
    exit 1
fi

# Execute the setup script
python3 -u semi_parametric_setup.py

echo "=== Done ==="
