#!/bin/sh

cd docs/
make clean
make html
cd ..
./clean_logs.sh
