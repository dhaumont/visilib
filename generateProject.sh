#! /bin/bash

mkdir debug
cd debug
cmake cmake -DCMAKE_BUILD_TYPE=Debug ..
make
cd ..

mkdir release
cd release
cmake cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

