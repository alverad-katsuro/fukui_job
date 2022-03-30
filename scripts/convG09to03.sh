#!/bin/bash
if [ "$#" -eq 0 ]; then echo -e "2010 by Lorenz Blum\nMakes gaussian09 output readable by molekel5.4\nUsage: $0 input output";exit;fi

#Gaussian 03: Otherwise it guesses it a Gaussian 94 outputfile, Density Matrix: To recognize the density matrix, Eigenvalues: To visualize the MOs, Atom AN: To see spectra
cat $1 | sed "s/Gaussian 09/Gaussian 03/" | sed "s/Eigenvalues -- /EIGENVALUES -- /" | sed "s/Density Matrix:/DENSITY MATRIX./" > $2
