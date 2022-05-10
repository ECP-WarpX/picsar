#!/usr/bin/env bash
#
# Copyright 2020 The PICSAR Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl, Luca Fedeli

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends \
    build-essential   \
    cmake             \
    libboost-dev      \
    libboost-math-dev \
    libboost-test-dev \
    g++ gfortran      \
    python3-dev

python3 -m pip install -U pip
python3 -m pip install "pybind11[global]"
