#!/bin/bash

for ((mN = 400 ; mN <= 1000 ; mN+=100)); do
  echo ./nu_cross.exe PDF4LHC21_mc $mN true ./
  ./nu_cross.exe PDF4LHC21_mc $mN true ./
done
