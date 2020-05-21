#! /bin/bash

# error handling
set -e -x

# platform definition
PLAT="manylinux1_x86_64"
# python version definition
PYVER="cp36-cp36m"  # possible options: (Python3.6)cp36-cp36m, (Python3.7)cp37-cp37m, (Python3.8)cp38-cp38m

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
PYBIN="/opt/python/${PYVER}/bin"
"${PYBIN}/pip" install -r /io/requirements.txt
"${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/


# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
  repair_wheel "$whl"
done