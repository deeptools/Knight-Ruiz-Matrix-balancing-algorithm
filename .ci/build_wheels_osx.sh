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
    $PYTHON setup.py bdist_wheel -d .
    delocate-wheel -w $WHEELHOUSE -v *.whl
    ls $WHEELHOUSE/krbalancing*.whl
    (
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
