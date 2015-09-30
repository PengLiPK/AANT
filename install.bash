#!/bin/bash

# Compile codes in src_anpro
cd src_anpro
rm *.o *.mod
make
cd ..

# Compile codes in src_groupvel
cd src_groupvel
rm *.o *.mod
make
cd ..

# Compile codes in src_tomo
cd src_tomo
rm *.o *.mod
make
cd ..
