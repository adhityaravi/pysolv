#!/bin/bash

# error handling
set -e -x

# platform definition
PLAT="manylinux1_x86_64"

# repair wheel
repair_wheel() {
    wheel="$1"
    if ! auditwheel show "$wheel"
    then
      echo "Skipping non-platform wheel $wheel"
    else
      auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

# Compile wheels
PYBIN="/opt/python/cp36-cp36m/bin"
"${PYBIN}/pip" install -r /io/requirements.txt
"${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/


# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
  repair_wheel "$whl"
done