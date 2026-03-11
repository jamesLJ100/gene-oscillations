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
import argparse

random.seed(42)
np.random.seed(42)
tf.random.set_seed(42)

MAX_LINEAR_DIMS = 5

def process_data_file(file_path, encoder_width, epochs, learning_rate, run_number):
    print(f"Processing file: {file_path}")
    print(f"Parameters: encoder_width={encoder_width}, epochs={epochs}, learning_rate={learning_rate}")
    
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
        epochs=epochs,
        rate=learning_rate,
        encoder_width=encoder_width,
        encoder_depth=len(encoder_width)
    )
    model.train(df.values, epochs=epochs, rate=learning_rate)
    
    elapsed = time.time() - start_time
    
    filepath, fullflname = os.path.split(file_path)
    fname, ext = os.path.splitext(fullflname)
    output_dir = os.path.join(filepath, 'cyclum')
    os.makedirs(output_dir, exist_ok=True)
    
    if run_number is not None:
        run_number_str = f"_r{run_number}"
        output_path = os.path.join(output_dir, f'{fname}{run_number_str}.h5')
    else:
        output_path = os.path.join(output_dir, f'{fname}.h5')
    

    mat2hdf(model.get_weight(), output_path)
    print(f"Model weights saved to: {output_path}")
    print(f"Runtime: {elapsed:.2f}s")
    
    return {"file": fname, "runtime_seconds": elapsed}

def main():
    parser = argparse.ArgumentParser(
        description='Process data files with Cyclum',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py input_dir --encoder_width 30 20 --epochs 500 --learning_rate 2e-4
  python script.py input_dir --encoder_width 40 30 20 --epochs 1000 --learning_rate 1e-3
        """
    )
    
    parser.add_argument('input_dir', type=str, help='Input directory containing H5 files')
    parser.add_argument('--encoder_width', nargs='+', type=int, default=[30, 20],
                        help='Encoder width layers (default: 30 20)')
    parser.add_argument('--epochs', type=int, default=500,
                        help='Number of training epochs (default: 500)')
    parser.add_argument('--learning_rate', type=float, default=2e-4,
                        help='Learning rate (default: 2e-4)')
    parser.add_argument('--run_number', type=int, default=None,
                        help='=Run Number')
    
    args = parser.parse_args()
    
    input_dir_abs = os.path.abspath(args.input_dir)
    print(f"Current working directory: {os.getcwd()}")
    print(f"Using parameters:")
    print(f"  Encoder width: {args.encoder_width}")
    print(f"  Epochs: {args.epochs}")
    print(f"  Learning rate: {args.learning_rate}")
    
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
        record = process_data_file(file_path, args.encoder_width, args.epochs, args.learning_rate, args.run_number)
        timing_records.append(record)
    
    timing_df = pd.DataFrame(timing_records)
    output_dir = os.path.join(input_dir_abs, 'cyclum')
    csv_path = os.path.join(output_dir, 'runtimes.csv')
    timing_df.to_csv(csv_path, index=False)
    print(f"\nRuntimes saved to: {csv_path}")

if __name__ == "__main__":
    main()
