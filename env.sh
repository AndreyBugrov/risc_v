#! /bin/bash

# $1 - device name

BASE= basename "$PWD"
if [ "$BASE" !=  "risc_v" ]; then
    git clone https://github.com/AndreyBugrov/risc_v
fi
git clone https://github.com/AndreyBugrov/OpenBLAS
mkdir open_blas
cd OpenBLAS
make
make install ../open_blas
sudo apt-get install libopenblas-dev
cd ..
mkdir bin
mkdir csv_results
g++ mat_opt/common.cpp mat_opt/main_test.cpp mat_opt/test.cpp mat_opt/multiplication.cpp -I open_blas/ -lopenblas -o bin/"$1"_test -O3 -fopenmp
