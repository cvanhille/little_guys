#!/bin/bash

#this script requires wget

STARTDIR=$(pwd)
cd "$( dirname "${BASH_SOURCE[0]}" )"
WDIR=$(pwd)

LAMMPSDIR=./lammps
MODFILESDIR=./LAMMPS_ModifiedFiles

echo "Downloading LAMMPS (3rd Nov 2022 stable release)..."

if [ ! -d $LAMMPSDIR ]; then
        if ! [ -x "$(command -v wget)" ]; then
                echo "You do not have wget installed on this computer, please download lammps 3Nov22 version manually and place it in a directory called 'lammps' at the same level as this file"
                echo "To do so:" 
                echo "  Copy this link into your browser: https://download.lammps.org/tars/lammps-3Nov2022.tar.gz"
                echo "  Extract this tar and copy the resulting file to the lammps directory"
        else
                wget -qO- https://download.lammps.org/tars/lammps-3Nov2022.tar.gz | tar xvz 
                mv lammps* $LAMMPSDIR
        fi
else
        echo "You already have a version of LAMMPS installed here. If you are unsure about which version it is delete the ./lammps directory and rerun this program"
        echo "To delete the current LAMMPS version just run 'rm -r ./lammps'"
fi

echo "Copying modified LAMMPS function files and compiling..."

if [ ! -f $LAMMPSDIR/src/lmp_serial ]; then
        cp -rf $MODFILESDIR/fix_bond_create.cpp $LAMMPSDIR/src/MC
        cp -rf $MODFILESDIR/fix_bond_create.h $LAMMPSDIR/src/MC
	cd $LAMMPSDIR/src
	make clean-all
	make yes-EXTRA-PAIR yes-MC yes-MISC yes-MOLECULE yes-REACTION yes-RIGID
	make serial
	cd "$WDIR"
else
        echo "Existing LAMMPS executable found! If you are unsure about which version it is delete the ./lammps directory and rerun this program"
        echo "To delete the current LAMMPS version just run 'rm -r ./lammps'"
fi

cd "${STARTDIR}"
echo "Done!"

