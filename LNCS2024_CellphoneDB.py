import scanpy as sc
import pandas as pd
import numpy as np
import glob
import warnings
import os
import sys

from cellphonedb.src.core.methods import cpdb_degs_analysis_method

pd.set_option('display.max_columns', 100)
# Define our base directory for the analysis
os.chdir('/cellphoneDB_analysis')

# Directories
cpdb_file_path = 'v5.0.0/cellphonedb.zip'
meta_file_path = 'metadata_hu_sLN_finalannotations.tsv'
counts_file_path = 'adata_rawcounts_cpdb.h5ad'
#microenvs_file_path = 'data/microenvironment.tsv'
degs_file_path = 'DEG_cellphoneDB.tsv' #taking from the DEGscran files
#active_tf_path = 'data/active_TFs.tsv'
out_path = 'results'

cpdb_results = cpdb_degs_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                            # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,                        # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
    counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
    #active_tfs_file_path = active_tf_path,                      # optional: defines cell types and their active TFs.
    #microenvs_file_path = microenvs_file_path,                  # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                                  # optional: whether to score interactions or not. 
    threshold = 0.1,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
    result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
    separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
    output_path = out_path,                                     # Path to save results
    output_suffix = None,                                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
    threads = 8
    )