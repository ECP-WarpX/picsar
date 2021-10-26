#!/usr/bin/env bash
#
# Copyright 2020 The PICSAR Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl, Luca Fedeli

set -eu -o pipefail

brew update
brew install boost
brew install pybind11
brew install libomp
#brew install open-mpi
