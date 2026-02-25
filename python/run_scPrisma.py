import pandas as pd
import os
import sys
import glob
import numpy as np
import scanpy as sc
import scPrisma
import h5py
import anndata as ad

ITERNUM = 100
REGULARISATION_STRENGTH = 0.1


def hdf2adata(file_path):
    with h5py.File(file_path, "r") as f:
        matrix = f["matrix"][:]
        rownames = f["rownames"][:].astype(str)
        colnames = f["colnames"][:].astype(str)
    
    adata = ad.AnnData(X=matrix)
    adata.obs.index = colnames
    adata.var.index = rownames
    
    return adata


def load_and_preprocess_data(file_path):
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
    
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")
    
    return adata


def apply_cyclic_reconstruction(adata):
    E_sga, E_rec_sga = scPrisma.algorithms_torch.reconstruction_cyclic_torch(
        adata.X, 
        iterNum=ITERNUM
    )
    
    sga_range = scPrisma.algorithms_torch.e_to_range(E_rec_sga)
    filtered_adata = adata[sga_range, :]
    
    print(f'E_rec_sga statistics - mean: {np.mean(E_rec_sga)}, std: {np.std(E_rec_sga)}')
    
    return filtered_adata, E_sga, E_rec_sga


def get_cyclic_gene_scores(adata):
    D = scPrisma.algorithms_torch.filter_cyclic_genes_torch(
        adata.X, 
        regu=REGULARISATION_STRENGTH,
        iterNum=ITERNUM
    )
    print(f'D matrix statistics - min={D.min()}, max={D.max()}, mean={D.mean()}')
    
    # Flip: D[i,i] close to 0 = cyclic, so 1 - D[i,i] close to 1 = cyclic
    cyclic_scores = 1 - np.diag(D)
    
    scores_df = pd.DataFrame({
        'symbol': adata.var.index,
        'score': cyclic_scores
    })
    
    return scores_df


def process_file(input_path, output_dir):
    print(f'Processing file: {input_path}')
    
    # Step 1: load and preprocess
    adata = load_and_preprocess_data(input_path)
    
    # Step 2: order cells along the cycle
    filtered_adata, _, _ = apply_cyclic_reconstruction(adata)
    
    # Step 3: score genes by cyclicity
    scores_df = get_cyclic_gene_scores(filtered_adata)
    
    fname = os.path.splitext(os.path.basename(input_path))[0]
    output_path = os.path.join(output_dir, f'{fname}_scPrisma.csv')
    
    scores_df.to_csv(output_path, index=False)
    print(f'Saved cyclic scores to: {output_path}')
    print(f'Top 10 most cyclic genes:\n{scores_df.head(10)}')


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_dir>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    input_dir_abs = os.path.abspath(input_dir)
    print(f"Current working directory: {os.getcwd()}")
    
    if os.path.exists(input_dir_abs):
        print(f"Directory contents: {os.listdir(input_dir_abs)}")
    
    output_dir = os.path.join(input_dir_abs, 'scPrisma')
    os.makedirs(output_dir, exist_ok=True)
    
    search_pattern = os.path.join(input_dir_abs, "*.h5")
    data_files = glob.glob(search_pattern)
    data_files = [f for f in data_files if not f.endswith("_sim.h5")]
    print(f"Found {len(data_files)} H5 files to process")
    
    if not data_files:
        print("No H5 files found.")
        sys.exit(1)
    
    for idx, file_path in enumerate(data_files, 1):
        print(f"\nProcessing file {idx}/{len(data_files)}: {os.path.basename(file_path)}")
        process_file(file_path, output_dir)


if __name__ == "__main__":
    main()
