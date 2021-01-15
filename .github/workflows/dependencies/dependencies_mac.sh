#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

brew uninstall python@2.7.17
brew untap local/openssl
brew untap local/python2
brew update
brew install boost
brew install cmake
brew install libomp
#brew install open-mpi
