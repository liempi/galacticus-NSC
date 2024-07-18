#!/bin/sh

export GALACTICUS_EXEC_PATH=/Users/liempi/Galacticus/DarkCore

export INSTALL_PATH1=/opt/local
export INSTALL_PATH2=/Users/liempi/Galacticus/packages


export FCCOMPILER=gfortran-12
export CCOMPILER=gcc-12
export CPPCOMPILER=g++-12

export GALACTICUS_FCFLAGS="-fintrinsic-modules-path $INSTALL_PATH1/finclude -fintrinsic-modules-path $INSTALL_PATH1/include -fintrinsic-modules-path $INSTALL_PATH1/include/gfortran -fintrinsic-modules-path $INSTALL_PATH1/lib/gfortran/modules -L$INSTALL_PATH1/lib -fintrinsic-modules-path $INSTALL_PATH2/finclude -fintrinsic-modules-path $INSTALL_PATH2/include -fintrinsic-modules-path $INSTALL_PATH2/include/gfortran -fintrinsic-modules-path $INSTALL_PATH2/lib/gfortran/modules -L$INSTALL_PATH2/lib"

export GALACTICUS_CFLAGS="-I/usr/local/include -I/opt/local/include -L/usr/local/finclude -L/opt/local/include"
export GALACTICUS_CPPFLAGS="-I/usr/local/include -I/opt/local/include -L/usr/local/finclude -L/opt/local/include"
make -j32  Galacticus.exe
