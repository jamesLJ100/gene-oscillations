import glob
import os
import sys
import time
import pandas as pd
import cyclum.tuning
import cyclum.models
from cyclum.hdfrw import mat2hdf, hdf2mat
import numpy as np
import sklearn as skl
import h5py
import random
import tensorflow as tf

random.seed(42)
np.random.seed(42)
tf.random.set_seed(42)

MAX_LINEAR_DIMS = 5
LEARNING_RATE = 2e-4
EPOCHS = 500

def process_data_file(file_path):
    print(f"Processing file: {file_path}")
    
    df = hdf2mat(file_path)
    print(f"Matrix shape: {df.shape}")  # N.B. cells should be rows
    print(f"First few row names: {df.index[:5].tolist()}")
    print(f"First few col names: {df.columns[:5].tolist()}")
    print(df.head())
    
    df = np.log2(df + 1)
    df = pd.DataFrame(
        skl.preprocessing.scale(df),
        index=df.index,
        columns=df.columns
    )
    
    start_time = time.time()
    
    model = cyclum.tuning.CyclumAutoTune(
        df.values,
        max_linear_dims=MAX_LINEAR_DIMS,
        early_stop=True,
        epochs=EPOCHS,
        rate=LEARNING_RATE,
        encoder_width=[30, 20]
    )
    model.train(df.values, epochs=EPOCHS, rate=LEARNING_RATE)
    
    elapsed = time.time() - start_time
    
    filepath, fullflname = os.path.split(file_path)
    fname, ext = os.path.splitext(fullflname)
    output_dir = os.path.join(filepath, 'cyclum')
    os.makedirs(output_dir, exist_ok=True)
    
    output_path = os.path.join(output_dir, f'{fname}_cyclum.h5')
    mat2hdf(model.get_weight(), output_path)
    print(f"Model weights saved to: {output_path}")
    print(f"Runtime: {elapsed:.2f}s")
    
    return {"file": fname, "runtime_seconds": elapsed}

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_dir>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    input_dir_abs = os.path.abspath(input_dir)
    print(f"Current working directory: {os.getcwd()}")
    
    if os.path.exists(input_dir_abs):
        print(f"Directory contents: {os.listdir(input_dir_abs)}")
    
    search_pattern = os.path.join(input_dir_abs, "*.h5")
    data_files = glob.glob(search_pattern)
    # exclude sim files
    data_files = [f for f in data_files if not f.endswith("_sim.h5")]
    print(f"Found {len(data_files)} H5 files to process")
    
    if not data_files:
        print("No H5 files found.")
        sys.exit(1)
    
    timing_records = []
    
    for idx, file_path in enumerate(data_files, 1):
        print(f"\nProcessing file {idx}/{len(data_files)}: {os.path.basename(file_path)}")
        record = process_data_file(file_path)
    
    timing_df = pd.DataFrame(timing_records)
    output_dir = os.path.join(input_dir_abs, 'cyclum')
    csv_path = os.path.join(output_dir, 'runtimes.csv')
    timing_df.to_csv(csv_path, index=False)
    print(f"\nRuntimes saved to: {csv_path}")

if __name__ == "__main__":
    main()
