#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Title: This module is generate report for PySeqRNA run.

:Created: May 15, 2022

:Author: naveen duhan

"""




from calendar import different_locale
from matplotlib.pyplot import title
import pandas as pd
import re
import os
import glob



def generate_report(outdir, combinations, infile, FOLD, FDR):

    title='Transcriptomics Analysis'
    with open(infile, 'r') as fp:
        lines = fp.readlines()

        for line in lines:
            if line.startswith("#"):
                title = line.replace("#", '')




    ### Check All directories 
    quality = os.path.join(outdir,'1_Quality')

    alignment = os.path.join(outdir,'2_Alignment')

    quantification = os.path.join(outdir,'3_Quantification')

    differential = os.path.join(outdir,'4_Differential_Expression')

    plots = os.path.join(outdir,'5_Visualization')

    annotation = os.path.join(outdir,'6_Functional_Annotation')

    ## Generate Quality Navbar
    quality_header='\n<li>\n<a href="#submenu1" data-bs-toggle="collapse" class="nav-link px-0 align-middle">\n<i class="fs-4 bi-speedometer2"></i> <span class="ms-1 d-none d-sm-inline">1. Quality</span> </a>\n<ul class="collapse show nav flex-column ms-1" id="submenu1" data-bs-parent="#menu">'
    quality_footer = '\n</ul>\n</li>'
    fastqc_1 = '\n<li class="w-100"><a href="#fastqc" class="nav-link pl-3"> <span class="d-none d-sm-inline">FastQC Results</span></a>\n</li>'
    fastqc_2 = '\n<li class="w-100"><a href="#trim_fastqc" class="nav-link pl-3"> <span class="d-none d-sm-inline">Trimmed FastQC Results</span></a>\n</li>'
    trim_galore = '\n<li class="w-100"><a href="#trim" class="nav-link pl-3"> <span class="d-none d-sm-inline">Trim Galore Results</span></a>\n</li>'
    trimmomatic = '\n<li class="w-100"><a href="#trim" class="nav-link pl-3"> <span class="d-none d-sm-inline">Trimmomatic Results</span></a>\n</li>'
    flexbar = '\n<li class="w-100"><a href="#trim" class="nav-link pl-3"> <span class="d-none d-sm-inline">Flexbar Results</span></a>\n</li>'

    if os.path.exists(os.path.join(quality, "fastqc_results")):

        quality_header +=fastqc_1

    if os.path.exists(os.path.join(quality, "trim_fastqc_results")):

        quality_header +=fastqc_2

    if os.path.exists(os.path.join(quality, "trim_galore_results")):

        quality_header +=trim_galore

    if os.path.exists(os.path.join(quality, "trimmomatic_results")):

        quality_header +=trimmomatic

    if os.path.exists(os.path.join(quality, "flexbar_results")):

        quality_header +=flexbar

    quality_header +=quality_footer

    ## Generate Alignment Navbar

    align_header='\n<li>\n<a href="#submenu2" data-bs-toggle="collapse" class="nav-link px-0 align-middle">\n<i class="fs-4 bi-speedometer2"></i> <span class="ms-1 d-none d-sm-inline">2. Alignment</span> </a>\n<ul class="collapse show nav flex-column ms-1" id="submenu2" data-bs-parent="#menu">'
    align_footer = '\n</ul>\n</li>'
    start_index = '\n<li class="w-100"><a href="#genome_index" class="nav-link pl-3"> <span class="d-none d-sm-inline">Genome Index</span></a>\n</li>'
    start_results = '\n<li class="w-100"><a href="#genome_alignment" class="nav-link pl-3"> <span class="d-none d-sm-inline">Aligned Reads</span></a>\n</li>'
    align_stats = '\n<li class="w-100"><a href="#astat" class="nav-link pl-3"> <span class="d-none d-sm-inline">Alignment Statistics</span></a>\n</li>'

    if os.path.exists(os.path.join(alignment, "star_index")):
        align_header +=start_index

    if os.path.exists(os.path.join(alignment, "hisat2_index")):
        align_header +=start_index

    if os.path.exists(os.path.join(alignment, "star_results")):
        align_header +=start_results

    if os.path.exists(os.path.join(alignment, "hisat2_results")):
        align_header +=start_results

    if os.path.exists(os.path.join(alignment, "alignment_statistics.xlsx")):
        align_header +=align_stats

    align_header += align_footer


    ## Generate Quantification Navbar
    quant_header= '\n<li>\n<a href="#submenu3" data-bs-toggle="collapse" class="nav-link px-0 align-middle">\n<i class="fs-4 bi-speedometer2"></i> <span class="ms-1 d-none d-sm-inline">3. Quantification</span> </a>\n<ul class="collapse show nav flex-column ms-1" id="submenu3" data-bs-parent="#menu">'
    quant_footer = '\n</ul>\n</li>'
    raw_counts = '\n<li class="w-100"><a href="#raw" class="nav-link pl-3"> <span class="d-none d-sm-inline">Raw Counts</span></a>\n</li>'
    rpkm_counts = '\n<li class="w-100"><a href="#rpkm" class="nav-link pl-3"> <span class="d-none d-sm-inline">RPKM Counts</span></a>\n</li>'
    tpm_counts = '\n<li class="w-100"><a href="#rpkm" class="nav-link pl-3"> <span class="d-none d-sm-inline">TPM Counts</span></a>\n</li>'
    cpm_counts = '\n<li class="w-100"><a href="#rpkm" class="nav-link pl-3"> <span class="d-none d-sm-inline">CPM Counts</span></a>\n</li>'
    median_counts = '\n<li class="w-100"><a href="#rpkm" class="nav-link pl-3"> <span class="d-none d-sm-inline">Median Ratio Counts</span></a>\n</li>'
    sample_cluster = '\n<li class="w-100"><a href="#cluster" class="nav-link pl-3"> <span class="d-none d-sm-inline">Sample Clustering</span></a>\n</li>'

    quant_header += raw_counts

    if os.path.exists(os.path.join(quantification, "RPKM_counts.xlsx")):
        quant_header += rpkm_counts
    if os.path.exists(os.path.join(quantification, "TPM_counts.xlsx")):
        quant_header +=tpm_counts
    if os.path.exists(os.path.join(quantification, "CPM_counts.xlsx")):
        quant_header +=cpm_counts
    if os.path.exists(os.path.join(quantification, "Median_ratio_counts.xlsx")):
        quant_header +=median_counts
    if os.path.exists(os.path.join(quantification, "Sample_cluster.png")):
        quant_header +=sample_cluster

    quant_header +=quant_footer

    ## Generate Differential Expression Navbar
    de_header= '\n<li>\n<a href="#submenu4" data-bs-toggle="collapse" class="nav-link px-0 align-middle">\n<i class="fs-4 bi-speedometer2"></i> <span class="ms-1 d-none d-sm-inline">4. Differential Expression</span> </a>\n<ul class="collapse show nav flex-column ms-1" id="submenu4" data-bs-parent="#menu">'
    de_footer = '\n</ul>\n</li>'

    all_gene = '\n<li class="w-100"><a href="#alldeg" class="nav-link pl-3"> <span class="d-none d-sm-inline">All Gene Expression</span></a>\n</li>'
    filt_gene = '\n<li class="w-100"><a href="#fdeg" class="nav-link pl-3"> <span class="d-none d-sm-inline">Filtered DEGs</span></a>\n</li>'
    diff_gene = '\n<li class="w-100"><a href="#diff_gene" class="nav-link pl-3"> <span class="d-none d-sm-inline">Differential Genes</span></a>\n</li>'

    de_header +=all_gene
    de_header +=filt_gene
    de_header +=diff_gene
    de_header +=de_footer


    ## Generate Plots Navbar
    plot_header= '\n<li>\n<a href="#submenu5" data-bs-toggle="collapse" class="nav-link px-0 align-middle">\n<i class="fs-4 bi-speedometer2"></i> <span class="ms-1 d-none d-sm-inline">5. Plots</span> </a>\n<ul class="collapse show nav flex-column ms-1" id="submenu5" data-bs-parent="#menu">'
    plot_footer = '\n</ul>\n</li>'

    maplot = '\n<li class="w-100"><a href="#ma" class="nav-link pl-3"> <span class="d-none d-sm-inline">MA Plots</span></a>\n</li>'
    volcano_plot = '\n<li class="w-100"><a href="#volcano" class="nav-link pl-3"> <span class="d-none d-sm-inline">Volcano Plots</span></a>\n</li>'
    venn_plot = '\n<li class="w-100"><a href="#venn" class="nav-link pl-3"> <span class="d-none d-sm-inline">VENN Plots</span></a>\n</li>'


    if os.path.exists(os.path.join(plots, 'MA_Plots')):
        plot_header += maplot

    if os.path.exists(os.path.join(plots, 'Volcano_Plots')):
        plot_header += volcano_plot

    if os.path.exists(os.path.join(plots, 'Venn_Plots')):
        plot_header += venn_plot

    plot_header +=plot_footer

    ## Generate Functional Annotation Navbar
    fa_header= '\n<li>\n<a href="#submenu6" data-bs-toggle="collapse" class="nav-link px-0 align-middle">\n<i class="fs-4 bi-speedometer2"></i> <span class="ms-1 d-none d-sm-inline">6. Functional Annotation</span> </a>\n<ul class="collapse show nav flex-column ms-1" id="submenu6" data-bs-parent="#menu">'
    fa_footer = '\n</ul>\n</li>'

    go = '\n<li class="w-100"><a href="#go" class="nav-link pl-3"> <span class="d-none d-sm-inline">Gene Ontology</span></a>\n</li>'
    kegg = '\n<li class="w-100"><a href="#kegg" class="nav-link pl-3"> <span class="d-none d-sm-inline">KEGG Pathway</span></a>\n</li>'


    if os.path.exists(os.path.join(annotation, 'Gene_Ontology')):
        fa_header += go

    if os.path.exists(os.path.join(annotation, 'KEGG_Pathway')):
        fa_header += kegg
    fa_header +=fa_footer


    copyright_header = '</ul>\n<hr>\n<div class="d-flex pb-4"><i class="fs-4 bi-speedometer2"></i> <span class="mx-1 d-none d-sm-inline">&copy; USU BioinfoCore 2022</span>  </div> </div> </div>'
    sample_input = pd.read_csv(infile, sep="\t", comment='#')
    infile = re.sub("class=\"dataframe ", "class=\"", pd.read_csv(infile, sep="\t", comment='#').to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))

    intro = f''' <div class="col-md-9 py-3 pe-5">
    <div class="row text-center py-2">
                        <h3>USU Bioinformatics Core Facility</h3>
                        <hr>

                    </div>
                    <div class="row text-center">
                        <h4>RNA Seq Analysis Report</h4>
                        <hr>

                    </div>
                    <div class="row justify-content-center">
    <div class="col-md-2">
    <h5 class="px-5">Overview</h5>
    <hr>
    </div>
                    </div>
                    <div  class="row" id="intro">

                    <p class="px-5">
                        <b>Experiment:</b>{title}<br> 
                        <b>Sequence type:</b> Single-end reads<br>
                        <b>Total Samples:</b> {sample_input.shape[0]} samples<br>
                        <b>Pairwise comparisons:</b>{combinations}<br>
                        <b>Results Directory:</b> <a href=".">{outdir}</a> <br>
                        <b>Reference:</b> Reference was downloaded from ENSEMBL: <a href="http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz" target="_blank">Mouse Reference Genome</a><br>
                        <b>Feature File:</b> Feature File was downloaded from ENSEMBL: <a href="http://ftp.ensembl.org/pub/release-107/gff3/mus_musculus/Mus_musculus.GRCm39.107.gff3.gz" target="_blank">Mouse Feature File</a><br>
                        <br>
                        Following is the input file used with Sample information:
                    
                    </p>
                
                    </div>
                    <div class="row justify-content-center">
                    
    <div class="col-md-6">
    {infile}
    </div>
                    </div>
                    <hr>
                    
                '''

    # Read Quality Results 
    final_quality_header =' <div class="row justify-content-center">\n<div class="col-md-2">\n<h5 class="px-5">Quality</h5>\n<hr>\n</div>\n</div>'

    quality1 = '''
    <div  class="row" id="fastqc">
    <h5 class="mb-3"> Quality check </h5>
    <p class='text-justify'>For the initial quality check of input fastq reads, <code>fastQC</code> was used.
    fastQC (<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a>)
    is a quality control analysis tool designed to spot potential problems in high throughput sequencing datasets.
    All the results are stored in <a href="1_Quality/fastqc_results/" target="_blank"> fastqc_results</a> in the 1_Quality/ subfolder of your results directory</p>
    </div>
    '''

    trim1 = '''<div  class="row" id="trim">
    <h5 class="mb-3">Read quality trimming and adapter removal </h5>
    <p class="text-justify">Read trimming and adapter removal was performed using <code>Trim Galore</code>.
    Trim Galore (<a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a>)
    is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation
    positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
    The trimmed fastq files are stored in <a href="1_Quality/trim_galore_results/" target="_blank">trim_galore_results</a> under the 1_Quality/ subdirectory. </p>
    </div>
    '''
    trim2 = '''<h5 class="mb-3">Read quality trimming and adapter removal </h5>
    <p class="text-justify">Read trimming and adapter removal was performed using <code>Trim Galore</code>.
    Trim Galore (<a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a>)
    is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation
    positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
    The trimmed fastq files are stored in <a href="1_Quality/trim_galore_results/" target="_blank">trim_galore_results</a> under the 1_Quality/ subdirectory. </p>
    '''
    trim3 = '''<h5 class="mb-3">Read quality trimming and adapter removal </h5>
    <p class="text-justify">Read trimming and adapter removal was performed using <code>Trim Galore</code>.
    Trim Galore (<a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a>)
    is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation
    positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
    The trimmed fastq files are stored in <a href="1_Quality/trim_galore_results/" target="_blank">trim_galore_results</a> under the 1_Quality/ subdirectory. </p>
    '''
    quality2 = ''' <div  class="row" id="trim_fastqc"><h5 class="mb-3" > Post trimming quality check </h5>
    <p class='text-justify'>For the quality check of trimmed fastq reads, <code>fastQC</code> was used.
    fastQC (<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a>)
    is a quality control analysis tool designed to spot potential problems in high throughput sequencing datasets.
    The results are stored in <a href="1_Quality/trim_fastqc_results/" target="_blank"> trim_fastqc_results</a> under the 1_Quality/ subdirectory.</p>
    </div>
    '''

    if os.path.exists(os.path.join(quality, "fastqc_results")):

        final_quality_header +=quality1

    if os.path.exists(os.path.join(quality, "trim_galore_results")):

        final_quality_header +=trim1

    if os.path.exists(os.path.join(quality, "trimmomatic_results")):

        final_quality_header +=trim2

    if os.path.exists(os.path.join(quality, "flexbar_results")):

        final_quality_header +=trim3

    if os.path.exists(os.path.join(quality, "trim_fastqc_results")):

        final_quality_header +=quality2

    final_quality_header +='<hr>'


    # Alignment results
    final_align_header =' <div class="row justify-content-center" >\n<div class="col-md-2">\n<h5 class="px-5" id="genome_index">Alignment</h5>\n<hr>\n</div>\n</div>'
    aligner1 = '''
    <h5 class="mb-3" id="genome_alignment">Reads alignment on reference genome</h5>
    <p clas="text-justify">To determine from where the reads are originated in the genome, the trimmed fastq reads were aligned with the reference genome using the <code>STAR</code> aligner. 
    STAR (<a href="https://github.com/alexdobin/STAR">https://github.com/alexdobin/STAR</a>) is an aligner designed to specifically address many of the challenges of RNA-Seq data mapping using a strategy to account for spliced alignments. 
    The alignments are stored in a <code>bam</code> format in <a href="2_Alignment/star_results/" target="_blank">star_results</a> under the 2_Alignments/ subfolder of your main output directory, Morrey_RNA-Seq. 
    Optionally, users can directly upload these <code>*.bam</code> files on any genomics viewer, like the <a href="https://software.broadinstitute.org/software/igv/home" target="_blank">Integrative Genomics Viewer</a> (IGV) for visualizing the RNA-Seq data. 
    IGV is a high-performance, easy-to-use, interactive tool for the visual exploration of genomics data.</p>
    '''
    aligner2 = '''
    <h5 class="mb-3">Reads alignment on reference genome</h5>
    <p clas="text-justify">To determine from where the reads are originated in the genome, the trimmed fastq reads were aligned with the reference genome using the <code>STAR</code> aligner. 
    STAR (<a href="https://github.com/alexdobin/STAR">https://github.com/alexdobin/STAR</a>) is an aligner designed to specifically address many of the challenges of RNA-Seq data mapping using a strategy to account for spliced alignments. 
    The alignments are stored in a <code>bam</code> format in <a href="2_Alignment/star_results/" target="_blank">star_results</a> under the 2_Alignments/ subfolder of your main output directory, Morrey_RNA-Seq. 
    Optionally, users can directly upload these <code>*.bam</code> files on any genomics viewer, like the <a href="https://software.broadinstitute.org/software/igv/home" target="_blank">Integrative Genomics Viewer</a> (IGV) for visualizing the RNA-Seq data. 
    IGV is a high-performance, easy-to-use, interactive tool for the visual exploration of genomics data.</p>
    '''


    stats = f'''<div class='row justify-content-center'> 
    <h5 class="mb-3" id="astat">Alignment statistics </h5>
    <p class="text-justify">The overall alignment statistics for each sample are recorded in <a href="2_Alignment/alignment_statistics.xlsx">alignment_statistics.xlsx</a> under the 2_Alignments/ subdirectory. This file contains 10 columns and all the samples. First column represents the Sample IDs; Second column represents the total number of input reads; Third column represents the total number of cleaned reads followed by the percentage values in the Fourth column; Fifth column represents the total number of reads aligned followed by the percent values in the Sixth column; Seventh column represents the uniquely mapped reads followed by the percentage in the Eighth column; and Ninth column represents the multi-mapped reads followed by the percentage values in the Tenth column.
    </p>'''

    if os.path.exists(os.path.join(alignment, "star_results")):
        final_align_header +=aligner1

    if os.path.exists(os.path.join(alignment, "hisat2_results")):
        final_align_header +=aligner2

    if os.path.exists(os.path.join(alignment, "alignment_statistics.xlsx")):
        final_align_header += stats
        final_align_header += '<div class="row justify-content-center"><div class="col-md-8" >'
        final_align_header += re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{alignment}/alignment_statistics.xlsx").to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_align_header += '</div></div>'
    final_align_header +='<hr>'


    # Quantification Results 
    final_quant_header = '<div class="row justify-content-center">\n<div class="col-md-2" id="raw">\n<h5 class="px-5">Quantification</h5>\n<hr>\n</div>\n</div>'

    fcount = f"""
    <h5 class="md-3"> Feature counts in the genome </h5>
    <p class="text-justify">This step counts reads for the given feature using <code>featureCounts</code>. FeatureCounts (<a href="http://subread.sourceforge.net/">http://subread.sourceforge.net/</a>) 
    is a program that counts how many reads map to genomic features, such as genes, exon, promoter and genomic bins. Raw aligner output is not sufficient for biological interpretation. Before read mapping results can be interpreted biologically,
    they must be summarized in terms of read coverage for genomic features of interest. Here, we perform feature counts for each sample. Results can be found in <a href="3_Quantification/Raw_Counts.xlsx" target="_blank"> Raw_counts.xlsx</a></p>
    """
    corr = """
    <h5 class="md-3" id="cluster"> Hierarchical cluster analysis of samples </h5>
    <p class="text-justify"> This step perform hierarchical cluster analysis of all the samples to check dissimilarity between samples. 
    Hierarchical clustering, also known as hierarchical clustering is a method for grouping similar objects into groups known as clusters. 
    The endpoint is a collection of clusters, each distinct from the others, and the objects within each cluster are broadly similar. 
    Raw counts can be used to perform hierarchical clustering leaving users with options to use RPKM values if needed. This graph can be downloaded here <a href="3_Quantification/Sample_cluster.png" target="_blank"> Sample_cluster</a> <br>
    <div class="row justify-content-center">
    <div class="col-md-6">
    <img src="3_Quantification/Sample_cluster.png"  height=600 width=auto></p>
    </div>
    <div>
    <hr>
    </div>
    </div>
    """

    if os.path.exists(os.path.join(quantification, "Raw_Counts.xlsx")):

        final_quant_header += fcount
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-8">'
        final_quant_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{quantification}/Raw_Counts.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_quant_header += '</div></div>'

    final_quant_header +=  f"""
    <h5 class="md-3" id="rpkm"> Normalized Counts </h5>
    <p class="text-justify">This step converts the raw counts into normalized counts. The numbers of mapped reads for each gene are proportional to the level of RNA expression, 
    which is both fascinating and uninteresting. Scaling raw count numbers to take into account the "uninteresting" elements is the process of normalization. In this manner,
     the expression levels within and/or between samples are more similar. 
    """
    if os.path.exists(os.path.join(quantification, "RPKM_counts.xlsx")):
        final_quant_header += 'We have normalized the raw read count using Reads Per Kilobase of transcript, per Million mapped reads (RPKM). The file is available at <a href="3_Quantification/RPKM_counts.xlsx">RPKM normalized counts</a> </p>'
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{quantification}/RPKM_counts.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "TPM_counts.xlsx")):
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{quantification}/TPM_counts.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "CPM_counts.xlsx")):
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{quantification}/CPM_counts.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "Median_ratio_counts.xlsx")):
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{quantification}/Median_ratio_counts.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "Sample_cluster.png")):
        final_quant_header += corr
       
   
    ## Differential Expression

    final_deg_header = '''<div class="row justify-content-center">\n<div class="col-md-4" id ="alldeg">\n<h5 class="px-5">Differential Expression</h5>\n<hr>\n</div>\n</div>'''
    final_deg_header += f'''<p>We have performed differential expression analysis of genes using DESeq2. There are 6 columns for each pairwise comparison, and  total {len(combinations)} pairwise comparisons. Scroll through this excel file (left to right). Here is how these 6 columns can be interpreted:<br>
    1. baseMean: It is a just the average of the normalized count values, dividing by size factors, taken over all the samples.<br> 
    2. logFC (log of Fold Change): logFC (fold change) generally refers to the ratio of average expression between two groups. For a particular gene, a log2 fold change of −1 for condition treated vs untreated means that the treatment induces a change in observed expression level of 2^−1 = 0.5 compared to the untreated condition. In simple terms, value 2 means that the expression has increased 4-fold, and so on. DESeq2 performs pair-wise tests for differential expression between two groups; log2FC=Log2(B)-Log2(A).
    <br>
    3. lfcSE: standard error of the log2FoldChange estimate.<br>

    4. stat = Wald statistic. By default, DESeq2 uses the Wald test to identify genes that are differentially expressed between two sample classes. The Wald statistic is the LFC divided by its standard error. This Wald statistic is used to calculate p-values (it is compared to a standard normal distribution). So, it's the ratio of LFC and SE which determines significance.
    
    <br>
    5. p-value: These are the Wald test p-values. <br>
    6. FDR: These are also called as the adjusted p-values (padj). By default, DESeq2 uses Benjamini-Hochberg method to adjust the p-values. </p>'''
    if os.path.exists(os.path.join(differential,'All_gene_expression.xlsx')):
        final_deg_header +=f'''\n Differential expression for all genes are presented in <a href="4_Differential_Expression/All_gene_expression.xlsx"> All gene expression</a>'''
        final_deg_header += '<div class="row justify-content-center"><div class="col-md-10 my-4">'
        final_deg_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{differential}/All_gene_expression.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_deg_header +='</div></div>'
    if os.path.exists(os.path.join(differential,'Filtered_DEGs.xlsx')):
        fd = pd.ExcelFile(os.path.join(differential,'Filtered_DEGs.xlsx'))
        fd_sheet = fd.sheet_names
        final_deg_header +=f'''\n Differential expression was filtered on user provided FOLD >= {FOLD} and FDR <= {FDR}. For example only one comparison is shown below. For each comparison there are different sheets in  <a href="4_Differential_Expression/Filtered_DEGs.xlsx"> Filtered differentially expressed genes</a> file. '''
        final_deg_header += '<div class="row justify-content-center"><div class="col-md-10 my-4" id="fdeg">'
        final_deg_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{differential}/Filtered_DEGs.xlsx",sheet_name=fd_sheet[0]).head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_deg_header +='</div></div>'
    if os.path.exists(os.path.join(differential,'Filtered_DEGs_summary.xlsx')):
        final_deg_header +=f'''\n Filtered DEGs summary is presented in <a href="4_Differential_Expression/Filtered_DEGs_summary.xlsx"> summary of differentially expressed genes</a>'''
        final_deg_header += '<div class="row justify-content-center"><div class="col-md-5 my-4">'
        final_deg_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{differential}/Filtered_DEGs_summary.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        
        final_deg_header +='</div></div>'
        final_deg_header +=f'''\n <p id="diff_gene">All filtered differentially expressed gene ID are present in <a href="./4_Differential_Expression/diff_genes/"> DEG gene IDs</a> </p>'''

    final_deg_header += f'''\n 
    Summary of filtered DEGs comparison wise are ploted:<br>
   <div class="row justify-content-center my-4">
    
    <div class="col-md-6 mt-2">
    <img src="./4_Differential_Expression/Filtered_DEG.png"  height=600></p>
    </div>
   </div>
    <hr>
    '''
  

  ## Plots 

    final_plots_header = f'''<div class="row justify-content-center">\n<div class="col-md-2">\n<h5 class="px-5">Visualization</h5>\n<hr>\n</div>\n</div>'''
    if os.path.exists(os.path.join(plots, 'Heatmap_top50.png')):
        final_plots_header += f'''<h5>Heatmap</h5> <br> Heatmaps are commonly used to visualize RNA-Seq results. They are useful for visualizing the expression of genes across the samples.
         The following heatmap was created basesd on top 50 ( 25 up, 25 down) differentially expressed genes.
         <div class="row justify-content-center my-4">
    
    <div class="col-md-6 mt-2">
    
    <img src="./5_Visualization/Heatmap_top50.png"  height=600></p>
    </div>
   </div>
         '''
    if os.path.exists(os.path.join(plots, 'MA_Plots')):

        maplots = glob.glob(f"{plots}/MA_Plots/*")
        final_plots_header += f'''<h5 id="ma">MA Plots</h5> <br> A 2-dimensional (2D) scatter plot called an MA plot is used to display gene expression datasets. 
        The MA plot uses the log of the mean of the normalized expression counts of the two conditions on the X-axis and the log of the fold change (M) 
        on the Y-axis to display and detect changes in gene expression from two distinct conditions (for example, normal vs. treated). In general, 
        log fold changes for genes with lower mean expression levels will be quite varied. Genes expressed similarly in both normal and treated samples 
        will group together around the M=0 value, i.e. genes expressed similarly across all treatments. Genes with considerable expression are shown by
        points away from the M=0 line. For instance, a gene is upregulated and downregulated if the point is above and below the M=0 line, respectively. Only one comparison plots is depicted below all other comparison MA plots are available 
        <a href="./5_Visualization/MA_Plots">MA Plots</a>
         <div class="row justify-content-center my-4">
    
        <div class="col-md-6 mt-2">
        
        <img src="{maplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>
            '''
    if os.path.exists(os.path.join(plots, 'Volcano_Plots')):

        vplots = glob.glob(f"{plots}/Volcano_Plots/*")
        final_plots_header += f'''<h5 id="volcano">Volcano Plots</h5> <br> A 2-dimensional (2D) scatter plot with a volcano-like form is called a volcano plot. The log fold change (X-axis) and negative log10 of
         the p value are used to visualize and identify statistically significant gene expression changes from two distinct experimental circumstances (e.g., normal vs. treated) (Y-axis). 
         The p value decreases when the Y-axis point is raised. Significant differences in gene expression between the two situations are shown by the larger dispersion of data points in the volcano plot. 
         It is simple to identify genes with substantial changes by visualizing the expression of hundreds of genes gathered from omics research (e.g., transcriptomics, genomics, and proteomics). 
         Only one comparison plots is depicted below all other comparison Volcano plots are available
        <a href="./5_Visualization/Volcano_Plots">Volcano Plots</a>
         <div class="row justify-content-center my-4">
    
        <div class="col-md-6 mt-2">
        
        <img src="{vplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>
            '''

    if os.path.exists(os.path.join(plots, 'Venn_Plots')):

        vennplots = glob.glob(f"{plots}/Venn_Plots/*")
        final_plots_header += f'''<h5 id="venn">Venn Plots</h5> <br> A Venn diagram is a diagram that shows all possible logical relations between a finite collection of different comparisons. 
        <a href="./5_Visualization/Venn_Plots">Venn Plots</a>
            <div class="row justify-content-center my-4">

        <div class="col-md-6 mt-2">
        
        <img src="{vennplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>
        <hr>
            '''

# Functional Annotation

    final_func_header = f'''<div class="row justify-content-center">\n<div class="col-md-4">\n<h5 class="px-5">Functional Annotation</h5>\n<hr>\n</div>'''
    if os.path.exists(os.path.join(annotation, 'Gene_Ontology')):
        gofiles = glob.glob(f"{annotation}/Gene_Ontology/GO_Files/*")
        goplots = glob.glob(f"{annotation}/Gene_Ontology/GO_Plots/*")
        final_func_header += f'''<h5 id="go">Gene Ontology</h5> <br> Gene Ontology enrichment analysis provides information on the function of genes. It is divided in three categories Biological process (BP),
        Molecular fucntion (MF), and Cellular component (CC). Gene Ontology results provides plots as well as files. GO enrichment results contains 10 columns. Gene Ontology files and plots can be found at can be found at <a href="./6_Functional_Annotation/Gene_Ontology/">Gene ontology </a> 
        <br>
        1. GO ID : Gene Ontology ID <br>
        2. GO Term : Gene Ontology term description. <br>
        3. Ontology : Ontology type <br>
        4. GeneRatio : Number of DEGs present in particular GO ID out of total DEGs <br>
        5. BgRatio : Number of genes represents the GO ID out of total genes in the genome. <br>
        6. Pvalues : P values <br>
        7. Counts : Number of DEGs <br>
        8. FDR : Flase discovery rate <br>
        9. Genes : DEGs IDs <br>
        10. logPvalues : Enrichment score calulated by <code> log<sub>10</sub>(1-Pvalues)</code>
        </div>''' 

        final_func_header += '<div class="row justify-content-center"><div class="col-md-10 my-4">'
        final_func_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(gofiles[0]).head(20).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_func_header +='</div></div>'
        final_func_header +=f'''<div class="row justify-content-center">
        
        <div class="col-md-8 mt-2">
        
        <img src="{goplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>'''
    if os.path.exists(os.path.join(annotation, 'KEGG_Pathway')):
        keggfiles = glob.glob(f"{annotation}/KEGG_Pathway/KEGG_Files/*")
        keggplots = glob.glob(f"{annotation}/KEGG_Pathway/KEGG_Plots/*")
        
        final_func_header += f'''<h5 id="kegg">KEGG Pathway</h5> <br> KEGG pathway is a database resource for understanding high level functions of genes. KEGG pathway results provides plots as well as files. KEGG enrichment results contains 9 columns.
        KEGG enrichment files and plots can be found at <a href="./6_Functional_Annotation/KEGG_Pathway/">KEGG results</a> 
        <br>
        1. Pathway ID : KEGG Pathway ID <br>
        2. Description : Pathway description. <br>
        3. Ontology : Ontology type <br>
        4. GeneRatio : Number of DEGs present in particular Pathway ID out of total DEGs <br>
        5. BgRatio : Number of genes represents the Pathway ID out of total genes in the genome. <br>
        6. Pvalues : P values <br>
        7. Counts : Number of DEGs <br>
        8. FDR : Flase discovery rate <br>
        9. Genes : DEGs IDs <br>
        10. logPvalues : Enrichment score calulated by <code> log<sub>10</sub>(1-Pvalues)</code>
   
        ''' 

        final_func_header += '<div class="row justify-content-center"><div class="col-md-10 my-4">'
        final_func_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(keggfiles[0]).head(20).to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_func_header +='</div></div></div>'
        final_func_header +=f'''<div class="row justify-content-center my-4">
    
        <div class="col-md-6 mt-2">
        
        <img src="{keggplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>
            '''

    header = '''<html>
    <head>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.0/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-gH2yIJqKdNHPEq0n4Mqa/HGKIhSkIHeL5AyhkYV8i59U5AR6csBvApHHNl/vI1Bx" crossorigin="anonymous">
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.0/dist/js/bootstrap.bundle.min.js" integrity="sha384-A3rJD856KowSb7dwlZdYEkO39Gagi7vIsF0jrRAoQmDKKtQBHUuLZ9AsSv4jD4Xa" crossorigin="anonymous"></script>
    <style>

    thead {
        background: #da4b6c;
    border-bottom: 1px solid #c9c4c4;
    color: #fff;
    }

    .table-responsived{
        display:block;
        width:100%;
        overflow-x: auto;
        max-height: 580px;
        border: none;
        padding: 5px 2px;
        
    }

    td{
        text-align:center!important;
    }

    .side{
        font-size:12px
    }

    </style>
    </head>
    <body>
    <div class="container-fluid">
    <div class="row ">
    <div class="col-md-3"></div>
    <div class="col-md-3 col-xl-2 px-0 bg-dark side fixed-top">
    <div class="d-flex flex-column align-items-center align-items-sm-start px-5 pt-5 text-white min-vh-100">
    <a href="" class="d-flex align-items-center pb-3 mb-md-0 me-md-auto text-white text-decoration-none">
    <span class="fs-5 d-none d-sm-inline">USU BioinfoCore</span>
    </a>
    <ul class="nav nav-pills flex-column mb-sm-auto mb-0 align-items-center align-items-sm-start" id="menu">
    <li class="nav-item">
    <a href="#intro" class="nav-link align-center px-0">
        <i class="fs-4 bi-house"></i> <span class="ms-1 d-none d-sm-inline">Overview</span>
    </a>
    </li>'''

    footer = """
    </div>
    </div>
        </body>
        </html>
    """



    with open(f"{outdir}/analysis_report.html", 'w') as fp:

        fp.write(header)
        fp.write(quality_header)
        fp.write(align_header)
        fp.write(quant_header)
        fp.write(de_header)
        fp.write(plot_header)
        fp.write(fa_header)
        fp.write(copyright_header)
        fp.write(intro)
        fp.write(final_quality_header)
        fp.write(final_align_header)
        fp.write(final_quant_header)
        fp.write(final_deg_header)
        fp.write(final_plots_header)
        fp.write(final_func_header)
        fp.write(footer)
        fp.close()



# generate_report("/home/naveen/Documents/Phd_work/example/pySeqRNA_results.3", ['A1-A6'], "/home/naveen/Documents/Phd_work/example/input_Sample.txt", 2, 0.05)