#!/bin/bash

echo "Fetching the data"
echo "There's a lot of it, so be patient (and get a decent "
echo "internet connection.."

if [ -f o3a_bbhpop_inj_info.hdf ]; then
    wget https://dcc-lho.ligo.org/public/0168/P2000217/002/o3a_bbhpop_inj_info.hdf
else
    echo "Found o3a_bbhpop_inj_info.hdf, skipping"
fi

if [ -d all_posterior_samples ]; then
    wget https://dcc.ligo.org/public/0169/P2000223/007/all_posterior_samples.tar
    tar xf all_posterior_samples.tar
    rm all_posterior_samples.tar
else
    echo "Found all_posterior_samples, skipping"
fi

