#!/bin/bash
set -e

# echo "=== Running semi_parametric_setup.ipynb ==="
# echo "Python: $(python3 --version)"
# echo "Working directory: $(pwd)"
# echo "Contents: $(ls)"
# echo "MG_FOLDER_PATH: $MG_FOLDER_PATH"
# echo ""

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
