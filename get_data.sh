#!/bin/bash

echo "Fetching the data"
echo "There's a lot of it, so be patient (and get a decent "
echo "internet connection.."

if [ -f o3a_bbhpop_inj_info.hdf ]; then
    echo "Found o3a_bbhpop_inj_info.hdf, skipping"
else
    wget https://dcc-lho.ligo.org/public/0168/P2000217/002/o3a_bbhpop_inj_info.hdf
fi

if [ -d all_posterior_samples ]; then
    echo "Found all_posterior_samples, skipping"
else
    wget https://dcc.ligo.org/public/0169/P2000223/007/all_posterior_samples.tar
    tar xf all_posterior_samples.tar
    rm all_posterior_samples.tar
fi

if [ -f GWTC-1_sample_release.tar.gz ]; then
    echo "Found GWTC-1_sample_release.tar.gz"
else
    wget https://dcc.ligo.org/public/0157/P1800370/005/GWTC-1_sample_release.tar.gz
    tar xf GWTC-1_sample_release.tar.gz -C all_posterior_samples/ --strip 1
fi
