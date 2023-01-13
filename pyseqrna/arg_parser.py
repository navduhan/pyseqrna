#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Title: Argument parser module for pySeqRNA

:Created: October 11, 2021

:Author: Naveen Duhan
"""

import argparse
from ast import parse
from textwrap import dedent
import pyseqrna.version
import sys

ver = pyseqrna.version.__version__
parser = argparse.ArgumentParser(description="""pyseqrna {} : a python based RNAseq data analysis package""".format(ver),
usage="""%(prog)s  input_file, samples_path, reference_genome, feature_file [options]""",
epilog="""Written by Naveen Duhan (naveen.duhan@usu.edu),
Kaundal Bioinformatics Lab, Utah State University,
Released under the terms of GNU General Public Licence v3""",    
formatter_class=argparse.RawTextHelpFormatter )

mandatory = parser.add_argument_group("Required positional arguments")
mandatory.add_argument(
    "input_file", type=str,
    help="""Tab-delimited file containing sample information""")

mandatory.add_argument(
    'samples_path', type=str,
    help="Directory for raw reads")

mandatory.add_argument(
    "reference_genome", type=str,
    help="Path to the reference genome file.")

mandatory.add_argument(
    "feature_file", type=str,
    help="Path to the GTF/GFF file ")


parser.add_argument("--version", action="version", version= 'pyseqrna (version {})'.format(ver), help= "Show version information and exit")

internal = parser.add_argument_group("Internal arguments")

internal.add_argument( "--source", dest="source", default='ENSEMBL', 
    help="Provide the source database of reference and feature file\n[default: ENSEMBL]")

internal.add_argument( "--taxid", dest='taxid', default=None, 
    help="Provide ncbi taxonomy id of the species")

internal.add_argument( "--species", dest='species', default=None, 
    help="Provide ncbi taxonomy id of the species")

internal.add_argument( "--organismType", dest='speciestype', default='plants', 
    help="Provide Organism class")

internal.add_argument( "--outdir",  default='pySeqRNA_results', 
    help="create output directory name to write results.\n[default: pySeqRNA_results] ")

internal.add_argument('--paired', dest='paired', action='store_true', default=False, help="Enable paired end functionality in pySeqRNA\n[default:False]")

internal.add_argument("--fastqc", action="store_true", default=False, dest="fastqc", help= "Enable initial quality check on raw reads with fastqc \n[default:True]")

internal.add_argument("--fastqcTrim", action="store_true", default=False, dest="fastqc2", help= "Enable quality check on trimmed reads with fastqc\n[default:False]")

internal.add_argument('--ribosomal',  dest='ribosomal', action='store_true', default=False, help="Enable removal of ribosomal RNA from reads\n[default:False]")

internal.add_argument('--rnadb',required='--ribosomal' in sys.argv, dest='rnadb', action='store_true', default=False, help="Enable removal of ribosomal RNA from reads\n[default:False]")

internal.add_argument('--multimappedGroups', dest='mmgg', action='store_true', default=False, help="Enable multimapped gene group quantification \n[default:False]")

internal.add_argument('--minReadCounts', dest='minreadcount', action='store_true', default=100, help="Minimum number of reads to consider per sample for MMG \n[default:100]")

internal.add_argument('--percentSample', dest='percentsample', action='store_true', default=0.5, help="Minimum number of reads to consider in percent of samples for MMG \n[default:0.5]")

internal.add_argument('--combination', dest='combination', type=str, nargs='+', default='all', 
    help="""Provide space separated combination of samples to 
be compared for differential expresion.
For example M1-A1 M1-V1 Z1-M1\n[Default:all]""")

internal.add_argument('--fdr', dest='fdr',default=0.05,type=float,
    help="False Discovery Rate threshold\n[default:0.05]" )

internal.add_argument('--fold', dest='fold', default=2,type=int,
    help="""FOLD change value for filtering DEGs. 
Remember pyseqrna performs log2 of the given value\n[default:2]""")

internal.add_argument('--noreplicate', dest='noreplicate',action='store_false', default=True ,
    help="""Execute Differential gene expression with no replicate""")

internal.add_argument('--normalizeCount', dest='normalizecount',default='RPKM', choices=['RPKM', 'TPM', 'CPM', 'medianRatiocount'], 
    help="Convert raw read counts to normalized counts\n[default:RPKM]")

internal.add_argument('--heatmap', dest='heatmap',action='store_true', default=False, 
    help="Create heatmap\n[default:False]")

internal.add_argument('--heatmapType', dest='heatmaptype', default='degs', choices=['counts', 'degs'],
    help="""Create heatmap based on selected choice \n[default: counts]""" )

internal.add_argument('--maPlot', dest='maplot',action='store_true', default=False, 
    help="Create MA plot\n[default:True]")

internal.add_argument('--volcanoPlot', dest='volcanoplot',action='store_true', default=False, 
    help="Create Volcano plot\n[default:True]")

internal.add_argument('--vennPlot',  dest='vennplot', default=False, action='store_true', 
    help="Enables venplot of differentially expressed genes.\n[default: False] ")

internal.add_argument('--vennCombinations',  dest='venncombination', type=str, nargs='+', default='random',
    help="""Provide space separated 2-4 combination of samples to 
be compared for differential expresion.For example M1-A1 M1-V1 Z1-M1\n[Default is to make random vennplot of 4 combinations]."""
)
internal.add_argument('--cluster', dest='cluster',default=True, action='store_true',
help="Cluster samples to find dissimilarities in data")

annotation= parser.add_argument_group("Functional annotation arguments")

annotation.add_argument('--geneOntology', action="store_true", default=False, dest="geneontology",
    help="""Enables gene ontology functional enrichment using BioMart""")

annotation.add_argument('--keggPathway',  action="store_true", default=False, dest="keggpathway",
    help="""Enables kegg pathway functional enrichment using KEGG""")

tools = parser.add_argument_group("External tool arguments")

tools.add_argument("--trimming", action='store', dest="trimming", type=str, default='trim_galore', choices=['flexbar', 'trimmomatic', 'trim_galore'], 
    help="Select a tool for quality based read trimming.\n[default: trim_galore]"  )
tools.add_argument('--aligner', dest='aligners', default='STAR', choices=['STAR', 'hisat2'], 
    help="Select a read alignment tool.\n[default: STAR]")
tools.add_argument('--quantTool', dest='quantification',default= 'featureCounts', choices=['featureCounts','Htseq'], 
    help= "Select a feature quantification tool.\n[default:featureCounts]")
tools.add_argument('--deTool', dest="detool", default='DESeq2', choices=['DESeq2', 'edgeR'],
        help="Select a tool for differential expression.\n[default:DESeq2]")

compute = parser.add_argument_group("Computation arguments")

compute.add_argument('--slurm', dest='slurm', action='store_true', default=False, help="Enable SLURM job scheduling on HPC\n[default:False]")

compute.add_argument('--threads', dest='threads', action='store', default= "80% of available CPU", help="Number of processors/threads to use\n[default:80%% of available CPU]")

compute.add_argument('--memory', dest='memory', action='store', default=16, help="Max memory to use (in GB)\n[default:16]")


