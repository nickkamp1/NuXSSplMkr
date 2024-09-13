#!/bin/bash

#./nu_cross_hnl.exe GRV98lo_patched 0300 true ./
#./nu_cross_hnl.exe GRV98lo_patched 0600 true ./
#./nu_cross_hnl.exe GRV98lo_patched 1000 true ./
#./nu_cross_hnl.exe GRV98lo_patched 0100 true ./

strings=(
    #0001500
    #0002000
    #0003000
    #0004000
    0005000
    0006000
    0007000
    0008000
    0009000
    0001000
    0015000
    0020000
    0030000
    0040000
    0050000
    0060000
    0070000
    0080000
    0090000
    0010000
    0150000
    0200000
    0300000
    0400000
    0500000
    0600000
    0700000
    0800000
    0900000
    0100000
    1500000
    2000000
    3000000
    4000000
    5000000
    6000000
    7000000
    8000000
    9000000
)
for hnl_mass in "${strings[@]}"; do
    mkdir "M_${hnl_mass}MeV"
    ./nu_cross_hnl.exe GRV98lo_patched "$hnl_mass" true ./
done

./nu_cross_hnl.exe GRV98lo_patched 1500 true ./
./nu_cross_hnl.exe GRV98lo_patched 2000 true ./
./nu_cross_hnl.exe GRV98lo_patched 3000 true ./
./nu_cross_hnl.exe GRV98lo_patched 4000 true ./
./nu_cross_hnl.exe GRV98lo_patched 5000 true ./
./nu_cross_hnl.exe GRV98lo_patched 6000 true ./
./nu_cross_hnl.exe GRV98lo_patched 7000 true ./
./nu_cross_hnl.exe GRV98lo_patched 8000 true ./
./nu_cross_hnl.exe GRV98lo_patched 9000 true ./
./nu_cross_hnl.exe GRV98lo_patched 10000 true ./
./nu_cross_hnl.exe GRV98lo_patched 15000 true ./
./nu_cross_hnl.exe GRV98lo_patched 20000 true ./
./nu_cross_hnl.exe GRV98lo_patched 30000 true ./
./nu_cross_hnl.exe GRV98lo_patched 40000 true ./
./nu_cross_hnl.exe GRV98lo_patched 50000 true ./
./nu_cross_hnl.exe GRV98lo_patched 60000 true ./
./nu_cross_hnl.exe GRV98lo_patched 70000 true ./
./nu_cross_hnl.exe GRV98lo_patched 80000 true ./
./nu_cross_hnl.exe GRV98lo_patched 90000 true ./
./nu_cross_hnl.exe GRV98lo_patched 100000 true ./
./nu_cross_hnl.exe GRV98lo_patched 150000 true ./
./nu_cross_hnl.exe GRV98lo_patched 200000 true ./
./nu_cross_hnl.exe GRV98lo_patched 300000 true ./
./nu_cross_hnl.exe GRV98lo_patched 400000 true ./
./nu_cross_hnl.exe GRV98lo_patched 500000 true ./
./nu_cross_hnl.exe GRV98lo_patched 600000 true ./
./nu_cross_hnl.exe GRV98lo_patched 700000 true ./
./nu_cross_hnl.exe GRV98lo_patched 800000 true ./
./nu_cross_hnl.exe GRV98lo_patched 900000 true ./
./nu_cross_hnl.exe GRV98lo_patched 1000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 1500000 true ./
./nu_cross_hnl.exe GRV98lo_patched 2000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 3000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 4000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 5000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 6000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 7000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 8000000 true ./
./nu_cross_hnl.exe GRV98lo_patched 9000000 true ./
