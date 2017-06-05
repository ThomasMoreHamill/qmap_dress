#!/bin/sh
#set -x

module unload craype-sandybridge
module use $(dirname $(dirname $PWD))/modulefiles
module load nbm-intel

export USELOCALLIB="Y"
export LIBLOCAL=$(dirname $(dirname $PWD))/lib

make
