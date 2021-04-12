#!/bin/bash
set -e 
set -x

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Place to store wheels.
WHEELHOUSE=${1-$HOME/wheelhouse}
echo "Path to store wheels : $WHEELHOUSE"
rm -rf $WHEELHOUSE && mkdir -p $WHEELHOUSE

PYDIR37=/opt/python/cp37-cp37m/
PYDIR38=/opt/python/cp38-cp38/
PYDIR39=/opt/python/cp39-cp39/

export PATH=$PYDIR39/bin:$PATH
$PYDIR39/bin/python3 -m pip install conan

for PYDIR in $PYDIR39 $PYDIR38 $PYDIR37; do
    PYTHON=$PYDIR/bin/python

    # dependencies
    $PYTHON -m pip install auditwheel pytest
    (
        cd ..
        rm -rf build  || echo "build does not exists"
	    PYLIB=$(ls -d $PYDIR/lib/python3.*)
	    PYINDIR=$(ls -d $PYDIR/include/python3.*)
        $PYTHON setup.py build_ext \
            --cmake-extra-args \
            "-DPython3_EXECUTABLE=$PYTHON -DPython3_INCLUDE_DIR=$PYINDIR -DPython3_LIBRARY=$PYLIB"
        $PYTHON setup.py bdist_wheel --skip-build -d .
        $PYTHON -m auditwheel repair *.whl -w $WHEELHOUSE

	# install and test it
        $PYTHON -m pip install *.whl
        $PYTHON -c 'import krbalancing; print(dir(krbalancing))'
        $PYTHON -c 'import krbalancing; print(krbalancing.__version__ )'

	# remove the old wheel
        rm -rf *.whl 
    )
done

PYTHON=$PYDIR38/bin/python
$PYTHON -m pip install twine

ls -lh $WHEELHOUSE/*.whl

# If successful, upload using twine.
if [ -n "$PYPI_PASSWORD" ]; then
    $PYTHON -m twine upload $WHEELHOUSE/krbalancing*.whl \
        --user __token__ \
        --password $PYPI_PASSWORD \
        --skip-existing
else
    echo "PYPI password is not set. Not uploading...."
fi
