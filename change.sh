#!/bin/bash

# Change the next line to your needs
# {a..b..c} == a, a+c, a+2c, ..., b
# c == 1 by default
for g2 in -{1..1}.{3..4}{0..9}; do

  make clean
  mv Constants.h Constants.old

  # Make sure that when using this script,
  # G2 == 0.f in Constants.h
  sed "s/G2 0.f/G2 ($g2)/g" Constants.old > Constants.h
  make

  # This line specifies how to call and where to move
  # the executable
  mv sim.x ../dist_t10_largeN/n100_t10_g${g2}.x
  mv Constants.old Constants.h

done

make clean
