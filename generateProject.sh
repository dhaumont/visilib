#! /bin/bash

export LEDA_DIR=$HOME/code/LEDA-7
mkdir -p debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
cd ..

mkdir -p release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

