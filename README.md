[![GitHub stars](https://img.shields.io/github/stars/Julien-Livet/Sym.svg)](https://github.com/Julien-Livet/Sym/stargazers)
[![GitHub issues](https://img.shields.io/github/issues/Julien-Livet/Sym.svg)](https://github.com/Julien-Livet/Sym/issues)
[![License](https://img.shields.io/github/license/Julien-Livet/Sym.svg)](LICENSE)
![CI](https://github.com/Julien-Livet/Sym/actions/workflows/build.yml/badge.svg)

# Sym

## General purpose

Sym is a symbolic tool to deal symbolic expressions like GiNaC in C++ or sympy in Python.

## Installation

```
sudo apt update
sudo apt install -y cmake make git libginac-dev libboost-all-dev
sudo apt install -y libgoogle-glog-dev libatlas-base-dev libeigen3-dev libsuitesparse-dev
cd ~
mkdir sym_ws
cd sym_ws
git clone https://github.com/Julien-Livet/Sym.git
cd Sym
mkdir build && cd build
cmake ..
make -j$(nproc)
ctest -V
```
