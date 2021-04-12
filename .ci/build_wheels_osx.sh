#!/bin/bash

set -e 
set -x

brew install cmake || echo "Failed to install cmake"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#export PATH=/usr/local/bin:$PATH

WHEELHOUSE=$HOME/wheelhouse
rm -rf $WHEELHOUSE && mkdir -p $WHEELHOUSE


# Always prefer brew version.
PYTHON=$(which python)

if [ ! -f $PYTHON ]; then
    echo "Not found $PYTHON"
    exit -1
fi

$PYTHON -m pip install setuptools --upgrade 
$PYTHON -m pip install wheel --upgrade
$PYTHON -m pip install delocate --upgrade
$PYTHON -m pip install twine  --upgrade
$PYTHON -m pip install pytest  --upgrade

( 
    mkdir -p _build_wheel && cd _build_wheel
    ls -ltr
    cmake ../.. -DPython3_EXECUTABLE=$PYTHON
    make -j4 
    make wheel
    ctest --output-on-failure
    delocate-wheel -w $WHEELHOUSE -v *.whl
    ls $WHEELHOUSE/krbalancing*.whl

    ## NOTE: I am contantly getting  the following error in venv.
    ## $ python -c 'import krbalancing; print(krbalancing.__version__ )'
    ## Fatal Python error: PyMUTEX_LOCK(gil->mutex) failed

    ## create a virtualenv and test this.
    ##VENV=$(pwd)/venv
    ##rm -rf $VENV
    (
        # $PYTHON -m venv $VENV
        # source $VENV/bin/activate
        # which python
        # now use venv pyhton.
        $PYTHON --version
        $PYTHON -m pip install $WHEELHOUSE/krbalancing*.whl
        $PYTHON -c 'import krbalancing; print(krbalancing.__version__ )'
        $PYTHON -m pip uninstall -y krbalancing
    )
)

if [ -n "$PYPI_PASSWORD" ]; then
    echo "Did you test the wheels?"
    $PYTHON -m twine upload \
        -u __token__ -p $PYPI_PASSWORD \
        --skip-existing \
        $WHEELHOUSE/krbalancing*.whl
fi
