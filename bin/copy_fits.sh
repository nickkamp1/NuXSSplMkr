#!/bin/bash

strings=(
   0000000
   0000100
#    0000200
#    0000300
#    0000400
#    0000500
#    0000600
#    0000700
#    0000800
#    0000900
#    0001000
#    0001100
#    0001200
#    0001300
#    0001400
#    0001500
#    0001600
#    0001700
#    0001800
#    0001900
#    0002000
#    0002500
#    0004000
#    0005000
#    0006000
#    0007000
#    0008000
#    0009000
#    0010000
#    0015000
#    0020000
    # 0030000
    # 0040000
    # 0050000
    # 0060000
    # 0070000
    # 0080000
    # 0090000
    # 0150000
    # 0200000
    # 0300000
    # 0400000
    # 0500000
    # 0600000
    # 0700000
    # 0800000
    # 0900000
    # 0100000
    # 1500000
    # 2000000
    # 3000000
    # 4000000
    # 5000000
    # 6000000
    # 7000000
    # 8000000
    # 9000000
)
outdir="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/SIREN/resources/CrossSections/DipoleHNLDISSplines/DipoleHNLDISSpline-v2.0/"
#outdir="/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/SIREN/resources/CrossSections/HNLDISSplines/HNLDISSplines-v2.0/"
for hnl_mass in "${strings[@]}"; do
    subdir="M_${hnl_mass}MeV"
    sigma_file="${subdir}/sigma-nu-N-em-GRV98lo_patched_central"
    sigmabar_file="${subdir}/sigma-nubar-N-em-GRV98lo_patched_central"
    dsigma_file="${subdir}/dsdxdy-nu-N-em-GRV98lo_patched_central"
    dsigmabar_file="${subdir}/dsdxdy-nubar-N-em-GRV98lo_patched_central"
    mkdir "${outdir}/${subdir}"
    if [ -f "${sigma_file}_v2.fits" ]; then
        cp "${sigma_file}_v2.fits" "${outdir}/${sigma_file}.fits"
    fi
    if [ -f "${sigmabar_file}_v2.fits" ]; then
        cp "${sigmabar_file}_v2.fits" "${outdir}/${sigmabar_file}.fits"
    fi
    if [ -f "${dsigma_file}_v2.fits" ]; then
        cp "${dsigma_file}_v2.fits" "${outdir}/${dsigma_file}.fits"
    fi
    if [ -f "${dsigmabar_file}_v2.fits" ]; then
        cp "${dsigmabar_file}_v2.fits" "${outdir}/${dsigmabar_file}.fits"
    fi
done
