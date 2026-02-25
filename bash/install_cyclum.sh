#!/bin/bash

conda create -n cyclum_env python=3.7 -y
conda activate cyclum_env
conda install -c conda-forge cudatoolkit=10.0 cudnn=7.6 tensorflow-gpu=1.15 pandas numpy scikit-learn scipy h5py -y
pip install git+https://github.com/KChen-lab/Cyclum.git
conda env export > cyclum_environment.yml
