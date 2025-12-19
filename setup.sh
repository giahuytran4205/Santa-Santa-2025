#!/bin/bash

apt install build-essential

mkdir build

mkdir output

cd build
cmake ..
cmake --build .