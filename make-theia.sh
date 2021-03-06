#!/bin/sh
#set -x

# ------------------------------------------------------------------------------ 
#
# From makefile...
#
#INCS =          -I. -I$(G2_INC4) $(NETCDF_INCLUDE) $(HDF4_INCLUDE) \
#                $(HDF5_INCLUDE)
#LIBS =          $(MDL_LIB4) $(W3NCO_LIB4) $(G2_LIB4) $(BACIO_LIB4) \
#                $(JASPER_LIB) $(PNG_LIB) $(Z_LIB) $(NETCDF_LDFLAGS_C) \
#                $(NETCDF_LDFLAGS_F) $(HDF4_LDFLAGS) $(HDF5_LDFLAGS) \
#                $(INTEL_OMP_LIB)
#
# ------------------------------------------------------------------------------ 

module load intel/15.3.187
module load netcdf
module load hdf5

export G2_INC4="-I/scratch3/NCEPDEV/nwprod/lib/incmod/g2_v2.5.0_4"
export NETCDF_INCLUDE="-I$NETCDF/include"
export HDF5_INCLUDE="-I$HDF5/include"
export HDF5_LDFLAGS="-L$HDF5/lib -lhdf5 -lhdf5_fortran"
export NETCDF_LDFLAGS_C="-L$NETCDF/lib -lnetcdf"
export NETCDF_LDFLAGS_F="-L$NETCDF/lib -lnetcdff"
export MDL_LIB4=/scratch3/NCEPDEV/mdl/Eric.Engle/oper/mos_shared/lib/libmdl_4.a
export W3NCO_LIB4=/scratch3/NCEPDEV/nwprod/lib/libw3nco_4.a
export G2_LIB4=/scratch3/NCEPDEV/nwprod/lib/libg2_v2.5.0_4.a
export BACIO_LIB4=/scratch3/NCEPDEV/nwprod/lib/libbacio_4.a
export JASPER_LIB=/scratch3/NCEPDEV/nwprod/lib/libjasper.a
export PNG_LIB=/scratch3/NCEPDEV/nwprod/lib/libpng.a
export Z_LIB=/scratch3/NCEPDEV/nwprod/lib/libz.a
export INTEL_OMP_LIB=/apps/intel/composer_xe_2015.3.187/compiler/lib/intel64/libiomp5.a
export COMP=ifort

make

if [ $? -eq 0 ]; then
   mv -v blend_precip_downscale_gammadress ../../exec/.
fi
