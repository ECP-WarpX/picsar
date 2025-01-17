#!/usr/bin/env bash
#
# Copyright 2023 The BLAST Community
#
# License: BSD-3-Clause-LBNL
# Authors: Luca Fedeli

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential   \
    cmake             \
    clang-15          \
    clang-tidy-15     \
    libblas-dev       \
    libc++-15-dev     \
    libboost-dev      \
    libboost-math-dev \
    libboost-test-dev \
    libomp-15-dev     \
    python3-dev

python3 -m pip install -U pip
python3 -m pip install "pybind11[global]"
