#!/bin/sh
set -x

# Clear environment
module purge

# Load nbm-intel which is do a lot
module use $(dirname $(dirname $PWD))/modulefiles
module load nbm-intel

# Set up some variables to in order to compile directly
# from this script  
export USELOCALLIB="Y"
export LIBLOCAL=$(dirname $(dirname $PWD))/lib
export COMP=ftn

make
