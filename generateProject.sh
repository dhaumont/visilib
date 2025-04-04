#! /bin/bash

mkdir debug
cd debug
cmake cmake -DCMAKE_BUILD_TYPE=Debug ..
make  2>&1 | tee ../build.log
cd ..

#mkdir release
#cd release
#cmake cmake -DCMAKE_BUILD_TYPE=Release ..
#make
#cd ..

