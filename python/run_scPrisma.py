import pandas as pd
import os
import sys
import glob
import numpy as np
import scanpy as sc
import scPrisma
import h5py
import anndata as ad
import argparse


def hdf2adata(file_path):
    with h5py.File(file_path, "r") as f:
        matrix = f["matrix"][:]
        rownames = f["rownames"][:].astype(str)
        colnames = f["colnames"][:].astype(str)
    
    adata = ad.AnnData(X=matrix)
    adata.obs.index = colnames
    adata.var.index = rownames
    
    return adata


def load_and_preprocess_data(file_path, n_top_genes=5000):
    adata = hdf2adata(file_path)
    
    print(f"Shape: {adata.shape}")
    print(f"obs (rows) first 5: {adata.obs.index[:5].tolist()}")
    print(f"var (cols) first 5: {adata.var.index[:5].tolist()}")
    
    sc.pp.filter_genes(adata, min_cells=1)
    print(f"After filter_genes: {adata.shape}")
    
    sc.pp.normalize_total(adata)
    print(f"After normalize - mean counts per cell: {adata.X.sum(axis=1).mean():.1f}")
    
    sc.pp.log1p(adata)
    print(f"After log1p - max value: {adata.X.max():.3f}")
    
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")
    
    return adata


def apply_cyclic_reconstruction(adata, iternum=100):
    E_sga, E_rec_sga = scPrisma.algorithms_torch.reconstruction_cyclic_torch(
        adata.X, 
        iterNum=iternum
    )
    
    sga_range = scPrisma.algorithms_torch.e_to_range(E_rec_sga)
    filtered_adata = adata[sga_range, :]
    
    print(f'E_rec_sga statistics - mean: {np.mean(E_rec_sga)}, std: {np.std(E_rec_sga)}')
    
    return filtered_adata, E_sga, E_rec_sga


def get_cyclic_gene_scores(adata, regularisation_strength=0.1, iternum=100):
    D = scPrisma.algorithms_torch.filter_cyclic_genes_torch(
        adata.X, 
        regu=regularisation_strength,
        iterNum=iternum
    )
    print(f'D matrix statistics - min={D.min()}, max={D.max()}, mean={D.mean()}')
    
    # Flip: D[i,i] close to 0 = cyclic, so 1 - D[i,i] close to 1 = cyclic
    cyclic_scores = 1 - np.diag(D)
    
    scores_df = pd.DataFrame({
        'symbol': adata.var.index,
        'score': cyclic_scores
    })
    
    return scores_df


def process_file(input_path, output_dir, args):
    print(f'Processing file: {input_path}')
    
    # Step 1: load and preprocess
    adata = load_and_preprocess_data(input_path, n_top_genes=args.n_top_genes)
    
    # Step 2: order cells along the cycle
    filtered_adata, _, _ = apply_cyclic_reconstruction(adata, iternum=args.iternum)
    
    # Step 3: score genes by cyclicity
    scores_df = get_cyclic_gene_scores(
        filtered_adata, 
        regularisation_strength=args.regularisation_strength,
        iternum=args.iternum
    )
    
    fname = os.path.splitext(os.path.basename(input_path))[0]
    
    # Add run number to filename if specified
    if args.run_number is not None:
        output_filename = f'{fname}_r{args.run_number}.csv'
    else:
        output_filename = f'{fname}.csv'
    
    output_path = os.path.join(output_dir, output_filename)
    
    scores_df.to_csv(output_path, index=False)
    print(f'Saved cyclic scores to: {output_path}')


def main():
    parser = argparse.ArgumentParser(
        description='Process data files with scPrisma',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py input_dir --iternum 200 --regularisation_strength 0.2
  python script.py input_dir --n_top_genes 3000 --run_number 1
        """
    )
    
    parser.add_argument('input_dir', type=str, help='Input directory containing H5 files')
    
    # scPrisma-specific parameters
    parser.add_argument('--iternum', type=int, default=100,
                        help='Number of iterations for scPrisma algorithms (default: 100)')
    parser.add_argument('--regularisation_strength', type=float, default=0.1,
                        help='Regularisation strength for cyclic gene filtering (default: 0.1)')
    parser.add_argument('--n_top_genes', type=int, default=5000,
                        help='Number of highly variable genes to select (default: 5000)')
    parser.add_argument('--run_number', type=int, default=None,
                        help='Run number to append to output filenames')
    
    args = parser.parse_args()
    
    # Validate input directory
    input_dir_abs = os.path.abspath(args.input_dir)
    print(f"Current working directory: {os.getcwd()}")
    print(f"Input directory: {input_dir_abs}")
    
    if not os.path.exists(input_dir_abs):
        print(f"Error: Directory '{input_dir_abs}' does not exist.")
        sys.exit(1)
    
    print(f"Directory contents: {os.listdir(input_dir_abs)}")
    
    # Create output directory
    output_dir = os.path.join(input_dir_abs, 'scPrisma')
    os.makedirs(output_dir, exist_ok=True)
    
    # Find H5 files
    search_pattern = os.path.join(input_dir_abs, "*.h5")
    data_files = glob.glob(search_pattern)
    data_files = [f for f in data_files if not f.endswith("_sim.h5")]
    
    print(f"\nFound {len(data_files)} H5 files to process")
    print(f"Parameters:")
    print(f"  - iternum: {args.iternum}")
    print(f"  - regularisation_strength: {args.regularisation_strength}")
    print(f"  - n_top_genes: {args.n_top_genes}")
    if args.run_number is not None:
        print(f"  - run_number: {args.run_number}")
    
    if not data_files:
        print("No H5 files found.")
        sys.exit(1)
    
    # Process each file
    for idx, file_path in enumerate(data_files, 1):
        print(f"\n{'='*60}")
        print(f"Processing file {idx}/{len(data_files)}: {os.path.basename(file_path)}")
        print(f"{'='*60}")
        process_file(file_path, output_dir, args)


if __name__ == "__main__":
    main()
