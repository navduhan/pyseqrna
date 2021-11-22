#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: naveen duhan

Run pyseqrna
"""

from pyseqrna import arg_parser
from pyseqrna import pyseqrna_utils as pu
from pyseqrna import quality_check as qc
from pyseqrna import quality_trimming as qt
from pyseqrna import  aligners as al
from pyseqrna import pyseqrna_stats as ps
from pyseqrna import quantification as quant
from pyseqrna import differential_expression as de
from pyseqrna import pyseqrna_plots as pp
from pyseqrna import ribosomal as ribo
from pyseqrna import multimapped_groups as mmg
import pandas as pd
import dill
import time
import os, sys
from waiting import wait

log = pu.PyseqrnaLogger(mode="w", log="analysis")

startTime= time.ctime()

log.info("Analysis started at %s", startTime)

def main():
    # Get all the options from the user
    options, unknownargs = arg_parser.parser.parse_known_args()   
    # Create directory for results
    outdir = pu.make_directory(options.outdir)
    # Read input samples from file
    input_data = pu.read_input_file(options.input_file,options.samples_path, paired = options.paired)

    samples = input_data['samples']
    targets = input_data['targets']
    
    if options.combination != 'all':

        combination = options.combination
    else:

        combination = input_data['combinations']
    
    # check if fastqc is turned on 

    log.info("Starting with read quality check")

    if options.fastqc:
        
        jobid, fastqc_results = qc.fastqcRun(sampleDict=samples,outDir=outdir, slurm=options.slurm, paired=options.paired)

        if options.slurm:
            for job in jobid:
                wait(lambda: pu.check_status(job), waiting_for="quality to finish")
                log.info(f"Quality check completed for job {job}")

            log.info("Read quality check completed succesfully")
        else:
        
            log.info("Read quality check completed succesfully")

    # Trimming

    log.info(f"Read trimming started with {options.trimming}")

    if options.trimming == 'flexbar':

        outtrim, jobidt = qt.flexbarRun(sampleDict=samples,pairedEND=options.paired, slurm=options.slurm, mem=options.memory,cpu=options.threads, outDir=outdir)
    
    elif options.trimming == 'trimmomatic':

        outtrim, jobidt = qt.trimmomaticRun(sampleDict=samples,slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=outdir,paired=options.paired)
        
    elif options.trimming == 'trim_galore':

        outtrim, jobidt = qt.trim_galoreRun(sampleDict=samples,slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=outdir,paired=options.paired)
        
    if options.slurm:
        for job in jobidt:
            wait(lambda: pu.check_status(job), waiting_for="trimming to finish")
            log.info(f"Trimming completed for job {job}")

        log.info("Read trimming completed succesfully")
    else:
    
        log.info("Read trimming completed succesfully")

    # Check if read quality check is activated after trimming
    
    if options.fastqc2:
        
        jobid, fastqc_results = qc.fastqcRun(sampleDict=samples,out="fastqc_results_after_trimming", outDir=outdir, slurm=options.slurm, paired=options.paired)

        if options.slurm:
            for job in jobid:
                wait(lambda: pu.check_status(job), waiting_for="quality to finish")
                log.info(f"Quality check for trimmed read completed for job {job}")
    
            log.info("Read quality check for trimmed reads completed succesfully")
        else:
        
            log.info("Read quality check for trimmed reads completed succesfully")

    
    # Removal of ribosomal RNA

    if options.ribosomal:

        log.info("Removing ribosomal RNA ")
        
        outtrim, jobs = ribo.sortmernaRun(outtrim,outdir, options.paired,cpu=options.threads, slurm=options.slurm)

        if options.slurm:

            wait(lambda: pu.check_status(jobs), waiting_for="ribosomal rna removal to finish")

            log.info("Ribsomal RNA removal completed succesfully")
        else:
        
            log.info("Ribsomal RNA removal completed succesfully")

    # Alignment Section
    
    log.info("Starting Alignment process")

    if options.aligners == 'STAR':

        aligner= al.STAR_Aligner(genome=options.reference_genome,  slurm=options.slurm,  outDir=outdir)

        log.info("Genome indexing started")

        jobida = aligner.build_index(mem=options.memory,cpu=options.threads)
    
        if options.slurm:

            wait(lambda: pu.check_status(jobida), waiting_for="genome indexing to finish")
            log.info("Genome index completed succesfully")
        else:
            log.info("Genome index completed succesfully")

        log.info("Checking genome index")


        if aligner.check_index():
            log.info("Genome index is valid")
        else:
            log.info("Genome index is not valid please create it again")
            sys.exit()
        log.info("Alignment started")

        outalign, jobalign = aligner.run_Alignment(outtrim, pairedEND=options.paired, mem=options.memory, cpu=options.threads)
    
    elif options.aligners == 'hisat2':

        aligner = al.hisat2_Aligner(genome=options.reference_genome,  slurm=options.slurm,  outDir=outdir)

        log.info("Genome indexing started")

        jobida= aligner.build_index(mem=options.memory,cpu=options.threads)
    
        if options.slurm:

            wait(lambda: pu.check_status(jobida), waiting_for="genome indexing to finish")
            log.info("Genome index completed succesfully")
        else:
            log.info("Genome index completed succesfully")

        log.info("Checking genome index")

        if aligner.check_index():
            log.info("Genome index is valid")
        else:
            log.info("Genome index is not valid please create it again")
            sys.exit()
        log.info("Alignment started")

        outalign, jobalign = aligner.run_Alignment(outtrim, pairedEND=options.paired, mem=options.memory, cpu=options.threads)

    if options.slurm:
        for job in jobalign:
            wait(lambda: pu.check_status(job), waiting_for="alignment to finish")
            log.info(f"Alignment completed for job {job}")
        log.info("Read alignment completed succesfully")
    else:

        log.info("Read alignment completed succesfully")
    
        # Alignment statistics section

    log.info("Calculating alignment statistics")
    
    try:
        align_stat = ps.align_stats.align_stats(samples,outtrim, outalign,pairedEND=options.paired)

        align_stat.to_excel(outdir+"/alignment_statistics.xlsx", index=False)

        log.info(f"alignment stats completed and written to {outdir}/alignment_statistics.xlsx")

    except Exception:

        log.error("Unable to calculate read alignment")
    
     
    log.info("Feature Count from aligned reads started")

    if options.quantification == 'featureCounts':

        fjob = quant.featureCount(gff=options.feature_file, bamDict=outalign,slurm=options.slurm, mem= options.memory, cpu=options.threads, outDir=outdir)

        if options.slurm:

            wait(lambda: pu.check_status(fjob), waiting_for="quantification to finish")

        log.info("feature counts completed and written in %s/Counts.txt",outdir)
    
    elif options.quantification == 'Htseq':

        fjob = quant.htseqCount(gff=options.feature_file, bamDict=outalign,slurm=options.slurm, mem= options.memory, cpu=options.threads, outDir=outdir)

        if options.slurm:

            wait(lambda: pu.check_status(fjob), waiting_for="quantification to finish")

        log.info("feature counts completed and written in %s/Counts.txt",outdir)
    
        # Multi mapped 
    
    if options.mmgg:
        log.info("Counting multimapped read grouops")
       
        mmg_count = mmg.countMMG(samples,outalign,options.feature_file+".bed")
        mmg_count.to_excel(os.path.join(outdir,"mmg_count.xlsx", index=False))
        log.info("counting of multimapped read group finished")


    
        
if __name__ == '__main__':
    main()
