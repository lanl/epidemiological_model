#!/bin/sh

cd docs/
make clean
make html
cd ..

open docs/_build/html/index.html
