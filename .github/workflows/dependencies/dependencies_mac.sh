#!/usr/bin/env bash
#
# Copyright 2020 The PICSAR Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl, Luca Fedeli

set -eu -o pipefail

# newer runner images do not ship GCC anymore
if brew list --formula gcc &>/dev/null; then
    brew unlink gcc
fi
brew update
brew install boost
brew install pybind11
brew install llvm
#brew install open-mpi
brew install libomp
brew link --force libomp
