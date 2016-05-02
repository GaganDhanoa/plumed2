#!/bin/bash
rm \#*
trjconv_mpi -f imd.trr -s imd.tpr -o 2new.gro -skip 100 -sep
python selectedAtoms.py
