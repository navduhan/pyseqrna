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

def _read_file(infile):
    if infile.endswith(".csv") or infile.endswith(".txt"):
        df = pd.read_csv(infile, comment=None, header=None, names=["line"])
    elif infile.endswith(".xlsx"):
        df = pd.read_excel(infile, comment=None, header=None, names=["line"])
        df = df.fillna(' ')
    else:
        raise ValueError("Unsupported file format:", infile)
    
    return df

def _nsamples(infile):
    if infile.endswith(".txt"):
        df = pd.read_csv(infile, sep="\s+", comment='#')
    elif infile.endswith(".csv"):
        df = pd.read_csv(infile, comment='#')
    elif infile.endswith(".xlsx"):
        df = pd.read_excel(infile, comment='#')
    else:
        raise ValueError("Unsupported file format:", infile)
    
    return df.shape[0]

def _process_infile(infile):
    """
    Reads a text, CSV, or Excel file, filters out comments, and formats it as an HTML table.

    Args:
        infile (str): Path to the input file.

    Returns:
        str: The formatted HTML table.
    """

    try:
        # Determine file type and read using appropriate method
        if infile.endswith(".txt"):
            df = pd.read_csv(infile, sep="\s+", comment='#')
        elif infile.endswith(".csv"):
            df = pd.read_csv(infile, comment='#')
        elif infile.endswith(".xlsx"):
            df = pd.read_excel(infile, comment='#')
        else:
            raise ValueError("Unsupported file format: {}".format(infile))

        # Remove unnecessary HTML class attribute by sub-setting (improved efficiency)
        formatted_df = df.to_html(index=False, classes='table', justify='center')
        return formatted_df

    except ValueError as e:
        print("Error:", e)
        return None

def generate_report(outdir, combinations, genome, gff, infile, FOLD, FDR):

    
    title='Transcriptomics Analysis'
    


    ### Check All directories 
    quality = os.path.join(outdir,'1_Quality')

    alignment = os.path.join(outdir,'2_Alignment')

    quantification = os.path.join(outdir,'3_Quantification')

    differential = os.path.join(outdir,'4_Differential_Expression')

    plots = os.path.join(outdir,'5_Visualization')

    annotation = os.path.join(outdir,'6_Functional_Annotation')

    ## Generate Quality Navbar
    quality_header= """
                        <li class="accordion">
                        <a href="#" data-toggle="collapse" data-target="#collapsequality" aria-expanded="false"
                            aria-controls="collapseOne" class="collapsible">
                            <span class="fa fa-trello mr-3" aria-hidden="true"></span>Quality
                        </a>
                        <div id="collapsequality" class="collapse" aria-labelledby="headingOne">
                            <div>
                                <ul>
"""
    quality_footer = """
                                </ul>
                            </div>
                        </div>
                    </li>

"""
    fastqc_1 = '<li><a href="#fastqc">FastQC Results</a></li>'
    fastqc_2 = '<li><a href="#fastqc_trim">Trimmed FastQC Results</a></li>'
    trim_galore = '<li><a href="#trim">Trim Galore Results</a></li>'
    trimmomatic = '<li><a href="#trim">Trimmomatic Results</a></li>'
    flexbar = '<li><a href="#trim">Flexbar Results</a></li>'

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
    align_header= """
                        <li class="accordion">
                        <a href="#" data-toggle="collapse" data-target="#collapseAlign" aria-expanded="false"
                            aria-controls="collapseOne" class="collapsible">
                            <span class="fa fa-trello mr-3" aria-hidden="true"></span>Alignment
                        </a>
                        <div id="collapseAlign" class="collapse" aria-labelledby="headingOne">
                            <div>
                                <ul>
    """
    align_footer = """
                                    </ul>
                                </div>
                            </div>
                        </li>

    """
    star_index = '<li><a href="#genome_index">STAR Index</a></li>'
    star_results = '<li><a href="#genome_alignment">STAR Alignment</a></li>'
    hisat2_index = '<li><a href="#genome_index">HISAT2 Index</a></li>'
    hisat2_results = '<li><a href="#genome_alignment">HISAT2 Alignment</a></li>'
    align_stats = '<li><a href="#align_stats">Alignment Statistics</a></li>'

    if os.path.exists(os.path.join(alignment, "star_index")):
        align_header +=star_index

    if os.path.exists(os.path.join(alignment, "hisat2_index")):
        align_header +=hisat2_index

    if os.path.exists(os.path.join(alignment, "star_results")):
        align_header +=star_results

    if os.path.exists(os.path.join(alignment, "hisat2_results")):
        align_header +=hisat2_results

    if os.path.exists(os.path.join(alignment, "alignment_statistics.xlsx")):
        align_header +=align_stats

    align_header += align_footer


    ## Generate Quantification Navbar

    quant_header= """
                        <li class="accordion">
                        <a href="#" data-toggle="collapse" data-target="#collapseQuant" aria-expanded="false"
                            aria-controls="collapseOne" class="collapsible">
                            <span class="fa fa-trello mr-3" aria-hidden="true"></span>Quantification
                        </a>
                        <div id="collapseQuant" class="collapse" aria-labelledby="headingOne">
                            <div>
                                <ul>
    """
    quant_footer = """
                                    </ul>
                                </div>
                            </div>
                        </li>

    """

    raw_counts = '\n<li><a href="#raw">Raw Counts</a>\n</li>'
    rpkm_counts = '\n<li><a href="#rpkm">RPKM Counts</a>\n</li>'
    tpm_counts = '\n<li><a href="#rpkm">TPM Counts</a>\n</li>'
    cpm_counts = '\n<li><a href="#rpkm">CPM Counts</a>\n</li>'
    median_counts = '\n<li><a href="#rpkm">Median Ratio Counts</a>\n</li>'
    sample_cluster = '\n<li><a href="#cluster">Sample Clustering</a>\n</li>'

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

    quant_header += quant_footer

    ## Generate Differential Expression Navbar
    de_header= """
                        <li class="accordion">
                        <a href="#" data-toggle="collapse" data-target="#collapseDe" aria-expanded="false"
                            aria-controls="collapseOne" class="collapsible">
                            <span class="fa fa-trello mr-3" aria-hidden="true"></span>Differential Expression
                        </a>
                        <div id="collapseDe" class="collapse" aria-labelledby="headingOne">
                            <div>
                                <ul>
    """
    de_footer = """
                                    </ul>
                                </div>
                            </div>
                        </li>

    """
    all_gene = '\n<li><a href="#alldeg">All Gene Expression</a>\n</li>'
    filt_gene = '\n<li><a href="#fdeg">Filtered Gene Expression</a>\n</li>'
    diff_sum = '\n<li><a href="#diff_sum">Differential Gene Summary</a>\n</li>'
    diff_gene = '\n<li><a href="#diff_gene">Differential Gene List</a>\n</li>'
  

    de_header += all_gene
    de_header += filt_gene
    de_header += diff_sum
    de_header += diff_gene
    de_header += de_footer


    ## Generate Plots Navbar

    plot_header= """
                        <li class="accordion">
                        <a href="#" data-toggle="collapse" data-target="#collapsePlot" aria-expanded="false"
                            aria-controls="collapseOne" class="collapsible">
                            <span class="fa fa-trello mr-3" aria-hidden="true"></span>Visualization
                        </a>
                        <div id="collapsePlot" class="collapse" aria-labelledby="headingOne">
                            <div>
                                <ul>
    """
    plot_footer = """
                                    </ul>
                                </div>
                            </div>
                        </li>

    """

    maplot = '\n<li><a href="#ma">MA Plots</a>\n</li>'
    volcano_plot = '\n<li><a href="#volcano">Volcano Plots</a>\n</li>'
    venn_plot = '\n<li><a href="#venn">VENN Plots</a>\n</li>'

    if os.path.exists(os.path.join(plots, 'MA_Plots')):
        plot_header += maplot

    if os.path.exists(os.path.join(plots, 'Volcano_Plots')):
        plot_header += volcano_plot

    if os.path.exists(os.path.join(plots, 'Venn_Plots')):
        plot_header += venn_plot

    plot_header +=plot_footer

    ## Generate Functional Annotation Navbar
    fa_header= """
                        <li class="accordion">
                        <a href="#func" data-toggle="collapse" data-target="#collapseFa" aria-expanded="false"
                            aria-controls="collapseOne" class="collapsible">
                            <span class="fa fa-trello mr-3" aria-hidden="true"></span>Functional Annotation
                        </a>
                        <div id="collapseFa" class="collapse" aria-labelledby="headingOne">
                            <div>
                                <ul>
    """
    fa_footer = """
                                    </ul>
                                </div>
                            </div>
                        </li>

    """

    go = '<li><a href="#go">Gene Ontology</a></li>'
    kegg = '<li><a href="#kegg">KEGG Pathway</a></li>'

    if os.path.exists(os.path.join(annotation, 'Gene_Ontology')):
        fa_header += go

    if os.path.exists(os.path.join(annotation, 'KEGG_Pathway')):
        fa_header += kegg
    fa_header +=fa_footer


   
    
    infiled = _process_infile(infile)
    intro = f''' 
    
        <div  class="row pinfo" id="intro">

            <p class="px-5 pinfo">
                <b>Experiment:</b>{title}<br> 
                <b>Sequence type:</b> Single-end reads<br>
                <b>Total Samples:</b> {_nsamples(infile)} samples<br>
                <b>Pairwise comparisons:</b>{combinations}<br>
                <b>Results Directory:</b> <a href=".">{outdir}</a> <br>
                <b>Reference:</b> Reference was downloaded from ENSEMBL: <a href="{genome}" target="_blank">Reference Genome</a><br>
                <b>Feature File:</b> Feature File was downloaded from ENSEMBL: <a href="{gff}" target="_blank">Mouse Feature File</a><br>
                <br>
                Following is the input file used with Sample information:
            </p>
                
        </div>
        <div class="row justify-content-center">
                    
            <div class="col-md-10">
                {infiled}
            </div>
        </div>
                    
    '''

    # Read Quality Results 
    final_quality_header =''' <div id="intro" class="row justify-content-center">
                <h2 class="mb-2">Quality Trimming and Adapter Removal</h2>

            </div><hr>'''

    quality1 = '''
    <div class="row justify-content-center px-5" id="fastqc">
    <h5 class="mb-3"> Raw Read Quality Check </h5>
    <p class='text-justify'>

    For the initial quality check of input fastq reads, <code>fastQC</code> was used.
    fastQC (<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a>)
    is a quality control analysis tool designed to spot potential problems in high throughput sequencing datasets.
    All the results are stored in <a href="1_Quality/fastqc_results/" target="_blank"> fastqc_results</a> in the <code>1_Quality/</code> subfolder of your results directory</p>
    </div>
    '''

    trim1 = '''<div  class="row justify-content-center px-5" id="trim">
    <h5 class="mb-3">Raw Read Quality Trimming and Adapter Removal </h5>
    <p class="text-justify">Read trimming and adapter removal was performed using <code>Trim Galore</code>.
    Trim Galore (<a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a>)
    is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation
    positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
    The trimmed fastq files are stored in <a href="1_Quality/trim_galore_results/" target="_blank">trim_galore_results</a> under the <code>1_Quality/</code> subdirectory. </p>
    </div>
    '''
    trim2 = '''<h5 class="mb-3">Raw Read Quality Trimming and Adapter Removal </h5>
    <p class="text-justify">Read trimming and adapter removal was performed using <code>Trim Galore</code>.
    Trim Galore (<a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a>)
    is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation
    positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
    The trimmed fastq files are stored in <a href="1_Quality/trim_galore_results/" target="_blank">trim_galore_results</a> under the <code>1_Quality/</code> subdirectory. </p>
    '''
    trim3 = '''<h5 class="mb-3">Raw Read Quality Trimming and Adapter Removal </h5>
    <p class="text-justify">Read trimming and adapter removal was performed using <code>Trim Galore</code>.
    Trim Galore (<a href="https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/</a>)
    is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation
    positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
    The trimmed fastq files are stored in <a href="1_Quality/trim_galore_results/" target="_blank">trim_galore_results</a> under the <code>1_Quality/</code> subdirectory. </p>
    '''
    quality2 = ''' <div  class="row" id="trim_fastqc"><h5 class="mb-3" > Post Trimming Reads Quality Check </h5>
    <p class='text-justify'>For the quality check of trimmed fastq reads, <code>fastQC</code> was used.
    fastQC (<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank">https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</a>)
    is a quality control analysis tool designed to spot potential problems in high throughput sequencing datasets.
    The results are stored in <a href="1_Quality/trim_fastqc_results/" target="_blank"> trim_fastqc_results</a> under the <code>1_Quality/</code> subdirectory.</p>
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
    final_align_header =''' <div id="intro" class="row justify-content-center">
                <h2 class="mb-2">Read Alignment</h2>

            </div><hr>'''
    
    aligner1 = '''<div class='row justify-content-center px-5'>
    <h5 class="mb-3" id="genome_alignment">Read Alignment on Reference Genome</h5>
    <p clas="text-justify">To determine from where the reads are originated in the genome, the trimmed fastq reads were aligned with the reference genome using the <code>STAR</code> aligner. 
    STAR (<a href="https://github.com/alexdobin/STAR">https://github.com/alexdobin/STAR</a>) is an aligner designed to specifically address many of the challenges of RNA-Seq data mapping using a strategy to account for spliced alignments. 
    The alignments are stored in a <code>bam</code> format in <a href="2_Alignment/star_results/" target="_blank">star_results</a> under the 2_Alignments/ subfolder of your main output directory, Morrey_RNA-Seq. 
    Optionally, users can directly upload these <code>*.bam</code> files on any genomics viewer, like the <a href="https://software.broadinstitute.org/software/igv/home" target="_blank">Integrative Genomics Viewer</a> (IGV) for visualizing the RNA-Seq data. 
    IGV is a high-performance, easy-to-use, interactive tool for the visual exploration of genomics data.</p></div>
    '''
    aligner2 = '''<div class='row justify-content-center px-5'>
    <h5 class="mb-3">Read Alignment on Reference Genome</h5>
    <p clas="text-justify">To determine from where the reads are originated in the genome, the trimmed fastq reads were aligned with the reference genome using the <code>HISAT2</code> aligner. 
    STAR (<a href="https://github.com/alexdobin/STAR">https://github.com/alexdobin/STAR</a>) is an aligner designed to specifically address many of the challenges of RNA-Seq data mapping using a strategy to account for spliced alignments. 
    The alignments are stored in a <code>bam</code> format in <a href="2_Alignment/star_results/" target="_blank">star_results</a> under the 2_Alignments/ subfolder of your main output directory, Morrey_RNA-Seq. 
    Optionally, users can directly upload these <code>*.bam</code> files on any genomics viewer, like the <a href="https://software.broadinstitute.org/software/igv/home" target="_blank">Integrative Genomics Viewer</a> (IGV) for visualizing the RNA-Seq data. 
    IGV is a high-performance, easy-to-use, interactive tool for the visual exploration of genomics data.</p></div>
    '''


    stats = f'''<div class='row justify-content-center px-5'> 
    <h5 class="mb-3" id="align_stats">Alignment Statistics </h5>
    <p class="text-justify">The overall alignment statistics for each sample are recorded in <a href="2_Alignment/alignment_statistics.xlsx">alignment_statistics.xlsx</a> under the 2_Alignments/ subdirectory. This file contains 10 columns and all the samples. First column represents the Sample IDs; Second column represents the total number of input reads; Third column represents the total number of cleaned reads followed by the percentage values in the Fourth column; Fifth column represents the total number of reads aligned followed by the percent values in the Sixth column; Seventh column represents the uniquely mapped reads followed by the percentage in the Eighth column; and Ninth column represents the multi-mapped reads followed by the percentage values in the Tenth column.
    </p></div>'''

    if os.path.exists(os.path.join(alignment, "star_results")):
        final_align_header += aligner1

    if os.path.exists(os.path.join(alignment, "hisat2_results")):
        final_align_header += aligner2

    if os.path.exists(os.path.join(alignment, "alignment_statistics.xlsx")):
        final_align_header += stats
        final_align_header += '<div class="row justify-content-center"><div class="col-md-12" >'
        aligndf =pd.read_excel(f"{alignment}/alignment_statistics.xlsx")
        align_table = aligndf.to_html(index=False, classes='table', justify='center')
        final_align_header +=f'{align_table}'
        # final_align_header += re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{alignment}/alignment_statistics.xlsx").to_html(index=False, classes='table table-responsived table-borderless align-middle', justify='center'))
        final_align_header += '</div></div>'
    final_align_header +='<hr>'


    # Quantification Results 
    final_quant_header = ''' <div id="intro" class="row justify-content-center">
                <h2 class="mb-2">Read Quantification</h2>

            </div><hr>'''

    fcount = f"""<div class='row justify-content-center px-5'>
    <h5 class="md-3"> Feature counts in the genome </h5>
    <p class="text-justify">This step counts reads for the given feature using <code>featureCounts</code>. FeatureCounts (<a href="http://subread.sourceforge.net/">http://subread.sourceforge.net/</a>) 
    is a program that counts how many reads map to genomic features, such as genes, exon, promoter and genomic bins. Raw aligner output is not sufficient for biological interpretation. Before read mapping results can be interpreted biologically,
    they must be summarized in terms of read coverage for genomic features of interest. Here, we perform feature counts for each sample. Results can be found in <a href="3_Quantification/Raw_Counts.xlsx" target="_blank"> Raw_counts.xlsx</a></p>
    </div>"""
    corr = """<div class='row justify-content-center px-5 mt-5'>
    <h5 class="md-3" id="cluster"> Hierarchical cluster analysis of samples </h5>
    <p class="text-justify"> This step perform hierarchical cluster analysis of all the samples to check dissimilarity between samples. 
    Hierarchical clustering, also known as hierarchical clustering is a method for grouping similar objects into groups known as clusters. 
    The endpoint is a collection of clusters, each distinct from the others, and the objects within each cluster are broadly similar. 
    Raw counts can be used to perform hierarchical clustering leaving users with options to use RPKM values if needed. This graph can be downloaded here <a href="3_Quantification/Sample_cluster.png" target="_blank"> Sample_cluster</a> <br>
   
    <div class="col-md-10">
    <img src="3_Quantification/Sample_cluster.png"  height=auto width==800></p>
    </div>
    <div>
    <hr>
    </div>
    </div>
    """

    if os.path.exists(os.path.join(quantification, "Raw_Counts.xlsx")):
        df = pd.read_excel(f"{quantification}/Raw_Counts.xlsx").head(50)
        final_quant_header += fcount

        final_quant_header += f'<div class="row justify-content-center"><div class="col-md-10">'

        html_table = df.to_html(classes='table', index=False, justify='center')

        # final_quant_header +=re.sub("class=\"dataframe ", "class=\"", pd.read_excel(f"{quantification}/Raw_Counts.xlsx").head(50).to_html(index=False, classes='table table-responsived table-borderless', justify='center'))
        final_quant_header += f'{html_table}'
        final_quant_header += '</div></div>'

    final_quant_header +=  f"""<div class='row justify-content-center px-5 mt-5'>
    <h5 class="md-3" id="rpkm"> Normalized Counts </h5>
    <p class="text-justify">This step converts the raw counts into normalized counts. The numbers of mapped reads for each gene are proportional to the level of RNA expression, 
    which is both fascinating and uninteresting. Scaling raw count numbers to take into account the "uninteresting" elements is the process of normalization. In this manner,
     the expression levels within and/or between samples are more similar <p></div>. 
    """
    if os.path.exists(os.path.join(quantification, "RPKM_counts.xlsx")):
        final_quant_header += 'We have normalized the raw read count using Reads Per Kilobase of transcript, per Million mapped reads (RPKM). The file is available at <a href="3_Quantification/RPKM_counts.xlsx">RPKM normalized counts</a> </p>'
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=f' {pd.read_excel(f"{quantification}/RPKM_counts.xlsx").head(50).to_html(index=False, classes="table", justify="center")}'
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "FPKM_counts.xlsx")):
        final_quant_header += 'We have normalized the raw read count using Reads Per Kilobase of transcript, per Million mapped reads (RPKM). The file is available at <a href="3_Quantification/RPKM_counts.xlsx">RPKM normalized counts</a> </p>'
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=f' {pd.read_excel(f"{quantification}/FPKM_counts.xlsx").head(50).to_html(index=False, classes="table", justify="center")}'
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "TPM_counts.xlsx")):
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header +=f' {pd.read_excel(f"{quantification}/TPM_counts.xlsx").head(50).to_html(index=False, classes="table", justify="center")}'
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "CPM_counts.xlsx")):
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header += f' {pd.read_excel(f"{quantification}/CPM_counts.xlsx").head(50).to_html(index=False, classes="table", justify="center")}'
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "Median_ratio_counts.xlsx")):
        final_quant_header += '<div class="row justify-content-center"><div class="col-md-10">'
        final_quant_header += f' {pd.read_excel(f"{quantification}/Median_ratio_counts.xlsx").head(50).to_html(index=False, classes="table", justify="center")}'
        final_quant_header += '</div></div>'
    if os.path.exists(os.path.join(quantification, "Sample_cluster.png")):
        final_quant_header += corr
       
   
    ## Differential Expression

    final_deg_header = ''' <div id="intro" class="row justify-content-center px-5">
                <h2 class="mb-2">Differential Expression</h2>

            </div><hr>'''
    final_deg_header += f''' <div id="intro" class="row justify-content-center px-5">
    <p>We have performed differential expression analysis of genes using DESeq2. There are 6 columns for each pairwise comparison, and  total {len(combinations)} pairwise comparisons. Scroll through this excel file (left to right). Here is how these 6 columns can be interpreted:<br></p>
    <code>
    1. baseMean: It is a just the average of the normalized count values, dividing by size factors, taken over all the samples.<br> 
    2. logFC (log of Fold Change): logFC (fold change) generally refers to the ratio of average expression between two groups. For a particular gene, a log2 fold change of −1 for condition treated vs untreated means that the treatment induces a change in observed expression level of 2^−1 = 0.5 compared to the untreated condition. In simple terms, value 2 means that the expression has increased 4-fold, and so on. DESeq2 performs pair-wise tests for differential expression between two groups; log2FC=Log2(B)-Log2(A).
    <br>
    3. lfcSE: standard error of the log2FoldChange estimate.<br>

    4. stat = Wald statistic. By default, DESeq2 uses the Wald test to identify genes that are differentially expressed between two sample classes. The Wald statistic is the LFC divided by its standard error. This Wald statistic is used to calculate p-values (it is compared to a standard normal distribution). So, it's the ratio of LFC and SE which determines significance.
    
    <br>
    5. p-value: These are the Wald test p-values. <br>
    6. FDR: These are also called as the adjusted p-values (padj). By default, DESeq2 uses Benjamini-Hochberg method to adjust the p-values.<br> </code></div>'''
    if os.path.exists(os.path.join(differential,'All_gene_expression.xlsx')):
        final_deg_header +=f'''\n<p class="px-5 mt-3">Differential expression for all genes are presented in <a href="4_Differential_Expression/All_gene_expression.xlsx"> All gene expression</a></p>'''
        final_deg_header += '<div class="row justify-content-center"><div class="col-md-12 my-4">'
        final_deg_header +=f'{pd.read_excel(f"{differential}/All_gene_expression.xlsx").head(50).to_html(index=False, classes="table table-responsive", justify="center")}'
        final_deg_header +='</div></div>'
    if os.path.exists(os.path.join(differential,'Filtered_DEGs.xlsx')):
        fd = pd.ExcelFile(os.path.join(differential,'Filtered_DEGs.xlsx'))
        fd_sheet = fd.sheet_names
        final_deg_header +=f'''\n Differential expression was filtered on user provided FOLD >= {FOLD} and FDR <= {FDR}. For example only one comparison is shown below. For each comparison there are different sheets in  <a href="4_Differential_Expression/Filtered_DEGs.xlsx"> Filtered differentially expressed genes</a> file. '''
        final_deg_header += '<div class="row justify-content-center px-5"><div class="col-md-9 my-4" id="fdeg">'
        final_deg_header += f'{pd.read_excel(f"{differential}/Filtered_DEGs.xlsx",sheet_name=fd_sheet[0]).head(50).to_html(index=False, classes="table table-responsive", justify="center")}'
        final_deg_header +='</div></div>'
    if os.path.exists(os.path.join(differential,'Filtered_DEGs_summary.xlsx')):
        final_deg_header +=f'''\n Filtered DEGs summary is presented in <a href="4_Differential_Expression/Filtered_DEGs_summary.xlsx"> summary of differentially expressed genes</a>'''
        final_deg_header += '<div class="row justify-content-center px-5"><div class="col-md-5 my-4">'
        final_deg_header += f'{pd.read_excel(f"{differential}/Filtered_DEGs_summary.xlsx").to_html(index=False, classes="table", justify="center")}'
        
        final_deg_header +='</div></div>'
        final_deg_header +=f'''\n <p id="diff_gene">All filtered differentially expressed gene ID are present in <a href="./4_Differential_Expression/diff_genes/"> DEG gene IDs</a> </p>'''

    final_deg_header += f'''\n 
    Summary of filtered DEGs comparison wise are ploted:<br>
   <div class="row justify-content-center px-5 my-4">
    
    <div class="col-md-6 mt-2">
    <img src="./4_Differential_Expression/Filtered_DEG.png"  height=600></p>
    </div>
   </div>
    <hr>
    '''
  

  ## Plots 

    final_plots_header = ''' <div id="intro" class="row justify-content-center px-5">
                <h2 class="mb-2">Visualization</h2>

            </div><hr>'''
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

    final_func_header = '''<div id="func" class="row justify-content-center px-5">
                <h2 class="mb-2">Functional Annotation</h2>

            </div><hr>'''
    
    if os.path.exists(os.path.join(annotation, 'Gene_Ontology')):
        gofiles = glob.glob(f"{annotation}/Gene_Ontology/GO_Files/*")
        goplots = glob.glob(f"{annotation}/Gene_Ontology/GO_Plots/*")
        final_func_header += f'''<div class="row justify-content-center px-5"><h5 id="go">Gene Ontology</h5> <br> <p>Gene Ontology enrichment analysis provides information on the function of genes. It is divided in three categories Biological process (BP), Molecular fucntion (MF), and Cellular component (CC). Gene Ontology results provides plots as well as files. GO enrichment results contains 10 columns. Gene Ontology files and plots can be found at <a href="./6_Functional_Annotation/Gene_Ontology/">Gene ontology </a> </p>
        <br>
        <code>
        1. GO ID : Gene Ontology ID <br>
        2. GO Term : Gene Ontology term description. <br>
        3. Ontology : Ontology type <br>
        4. GeneRatio : Number of DEGs present in particular GO ID out of total DEGs <br>
        5. BgRatio : Number of genes represents the GO ID out of total genes in the genome. <br>
        6. Pvalues : P values <br>
        7. Counts : Number of DEGs <br>
        8. FDR : Flase discovery rate <br>
        9. Genes : DEGs IDs <br>
        10. logPvalues : Enrichment score calulated by  log<sub>10</sub>(1-Pvalues)</code>
        </div>''' 

        final_func_header += '<div class="row justify-content-center"><div class="col-md-10 my-4">'
        final_func_header +=f'{pd.read_excel(gofiles[0]).head(20).to_html(index=False, classes="table table-responsive", justify="center")}'
        final_func_header +='</div></div>'
        final_func_header +=f'''<div class="row justify-content-center">
        
        <div class="col-md-8 mt-2">
        
        <img src="{goplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>'''
    if os.path.exists(os.path.join(annotation, 'KEGG_Pathway')):
        keggfiles = glob.glob(f"{annotation}/KEGG_Pathway/KEGG_Files/*")
        keggplots = glob.glob(f"{annotation}/KEGG_Pathway/KEGG_Plots/*")
        
        final_func_header += f'''<div class="row justify-content-center px-5"><h5 id="kegg">KEGG Pathway</h5> <br> <p> KEGG pathway is a database resource for understanding high level functions of genes. KEGG pathway results provides plots as well as files. KEGG enrichment results contains 9 columns.
        KEGG enrichment files and plots can be found at <a href="./6_Functional_Annotation/KEGG_Pathway/">KEGG results</a> </p>
        <br>
        <code>
        1. Pathway ID : KEGG Pathway ID <br>
        2. Description : Pathway description. <br>
        3. Ontology : Ontology type <br>
        4. GeneRatio : Number of DEGs present in particular Pathway ID out of total DEGs <br>
        5. BgRatio : Number of genes represents the Pathway ID out of total genes in the genome. <br>
        6. Pvalues : P values <br>
        7. Counts : Number of DEGs <br>
        8. FDR : Flase discovery rate <br>
        9. Genes : DEGs IDs <br>
        10. logPvalues : Enrichment score calulated by  log<sub>10</sub>(1-Pvalues)</code>
   
        </div>''' 

        final_func_header += '<div class="row justify-content-center"><div class="col-md-10 my-4">'
        final_func_header += f'{pd.read_excel(keggfiles[0]).head(20).to_html(index=False, classes="table table-responsive", justify="center")}'
        final_func_header +='</div></div>'
        final_func_header +=f'''<div class="row justify-content-center my-4">
    
        <div class="col-md-6 mt-2">
        
        <img src="{keggplots[0].split(outdir+"/")[1]}"  height=600></p>
        </div>
        </div>
            '''

    html = f"""
            <!doctype html>
        <html lang="en">

       <head>
    <title>pySeqRNA</title>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css">
    <link rel="stylesheet" type="text/css"  href="https://cdn.jsdelivr.net/gh/navduhan/pyseqranhtml@main/style.css" crossorigin="anonymous">
    <link rel = "stylesheet" type="text/css" href= "https://cdn.jsdelivr.net/gh/navduhan/pyseqranhtml@main/custom.css">
     <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
   


  
    <meta name="robots" content="noindex, follow">

     </head>

        <body>
            <div class="wrapper d-flex align-items-stretch">
                <nav id="sidebar">
                    <div class="custom-menu">
                        <button type="button" id="sidebarCollapse" class="btn btn-primary">
                            <i class="fa fa-bars"></i>
                            <span class="sr-only">Toggle Menu</span>
                        </button>
                    </div>
            <div class="p-4">
                        <h1><a href="https://kaabil.net/pyseqrna" target="_blank" class="logo">pySeqRNA <span>version 0.2</span></a></h1>
                        <ul class="list-unstyled components mb-5">
                            <li class="active">
                                <a href="#intro"><span class="fa fa-home mr-3"></span>Introduction</a>
                            </li>
                            {quality_header}
                            {align_header}
                            {quant_header}
                            {de_header}
                            {plot_header}
                            {fa_header}
                    </ul>

                <div class="footer">
                    <p>
                        Copyright &copy;
                        <script>document.write(new Date().getFullYear());</script> pySeqRNA <i class="icon-heart"
                            aria-hidden="true"></i>
                    </p>
                </div>
            </div>
        </nav>
        <div id="content" class="p-4 p-md-5 pt-5" style="overflow-y: auto; max-height: calc(100vh - 100px);">
                    <div class="row justify-content-center">
                        <h2 class="mb-2">pySeqRNA Analysis Report </h2>

                    </div>
                    <hr>
                    <div id="intro" class="row justify-content-center">
                        <h2 class="mb-2">Analysis Information</h2>

                    </div>
                    <hr>
                <!-- Information seaction -->
                {intro}
                <hr>
<!-- quality section -->
{final_quality_header}

{final_align_header}
{final_quant_header}
{final_deg_header}
{final_plots_header}
{final_func_header}



        </div>
        </div>
        
    <script src="https://cdn.jsdelivr.net/gh/navduhan/pyseqranhtml@main/jquery.min.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/navduhan/pyseqranhtml@main/popper.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/navduhan/pyseqranhtml@main/bootstrap.min.js"></script>
    <script src="https://cdn.jsdelivr.net/gh/navduhan/pyseqranhtml@main/main.js"></script>

    <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
<script>
$(document).ready(function() {{
    $('.table').DataTable();
}});
</script>
<footer>
<div class="row justify-content-center">
Copyright © 2024 pySeqRNA
</div>

</footer>
        </body>

        </html>
                """


    with open(f"{outdir}/analysis_report.html", 'w') as fp:

        fp.write(html)

        fp.close()



# generate_report("/home/naveen/workspace/pySeqRNA_results.1", ['A1-M6'], "genome", "gff" ,"/home/naveen/workspace/input.txt", 2, 0.05)

