#!/bin/bash
python setup.py build
#python setup.py install
cp build/lib.*/_fastscoreCC.pyd _fastscoreCC.pyd
python test-fastscoreCC.py