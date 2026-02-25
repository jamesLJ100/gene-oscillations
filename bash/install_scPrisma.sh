#!/bin/bash
conda create -n scPrisma_env python=3.10 -y
conda activate scPrisma_env

git clone https://github.com/nitzanlab/scPrisma.git
cd scPrisma
pip install ".[gpu]"

pip install \
    pandas \
    numpy \
    scanpy \
    anndata \
    scipy \
    matplotlib \
    scikit-learn \
    h5py