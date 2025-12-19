#!/bin/bash

apt install build-essential -y
apt install cmake -y

mkdir build

mkdir output

cd build
cmake ..
cmake --build .