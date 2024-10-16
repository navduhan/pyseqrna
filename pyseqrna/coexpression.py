#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
:Title: This script contains clust function for co-expression of genes.

:Created: Feb 28, 2024

:Author: Naveen Duhan
"""

import os
import shutil
import glob
import sys
import pandas as pd
import subprocess
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu

# Initialize logger

log = PyseqrnaLogger(mode="a", log='coexpression')

def run_Clusts(countFile=None, targets=None, outdir="."):
    """_summary_

    Args:
        df (_type_): _description_
    """

    if countFile.endswith(".xlsx"):
          counts = pd.read_excel(countFile)
          counts.to_csv(f"{os.getcwd()}/raw_counts.txt", index=False, sep="\t")

          countFile = f"{os.getcwd()}/raw_counts.txt"
    else:
          countFile = countFile

    df = targets.groupby('sample')[0].apply(lambda x: ' '.join(x)).reset_index()

    df.to_csv(f"{os.getcwd()}/replicate_structure_file_clusts.txt", sep="\t", index=False)
    
    clustcmd = f"clust -n {countFile} -r {os.getcwd()}/replicate_structure_file_clusts.txt -o {outdir}"
    
    try:
        with open(os.path.join(outdir, "clust.out"), 'w+') as fout:
            with open(os.path.join(outdir, "clust.err"), 'w+') as ferr:

                job = subprocess.call(
                    clustcmd, shell=True, stdout=fout, stderr=ferr)

                log.info(
                    "Job successfully completed for gene co-expression")

    except Exception:

                    log.info("Job submition failed for gene co-expression ")

    os.remove(f"{os.getcwd()}/raw_counts.txt")
    os.remove(f"{os.getcwd()}/replicate_structure_file_clusts.txt")
    
    return
