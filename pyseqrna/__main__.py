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
from pyseqrna import gene_ontology as go
import pyseqrna.version
import pandas as pd
import numpy as np
import dill
import time
import os, sys
from waiting import wait

log = pu.PyseqrnaLogger(mode="w", log="analysis")



def main():
    
    # Get all the options from the user
    options, unknownargs = arg_parser.parser.parse_known_args()  

    if options.threads =='80% of available CPU' :
    
        options.threads = pu.get_cpu()
        
    startTime= time.ctime()


    log.info("Analysis started at %s", startTime)
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
        
        jobid, fastqc_results = qc.fastqcRun(sampleDict=samples,outDir=outdir, slurm=options.slurm, mem=options.memory, cpu=options.threads, paired=options.paired)

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

        outtrim, jobidt = qt.flexbarRun(sampleDict=samples,paired=options.paired, slurm=options.slurm, mem=options.memory,cpu=options.threads, outDir=outdir)
    
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
        
        jobid, fastqc_results = qc.fastqcRun(sampleDict=samples,out="fastqc_results_after_trimming", outDir=outdir, slurm=options.slurm, mem=options.memory, cpu=options.threads, paired=options.paired)

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
        
        outtrim, jobs = ribo.sortmernaRun(outtrim, outdir, rnaDatabases=options.rnadb, pairedEND= options.paired,cpu=options.threads, slurm=options.slurm)

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
        align_stat = ps.align_stats(samples,outtrim, outalign,pairedEND=options.paired)

        align_stat.to_excel(outdir+"/alignment_statistics.xlsx", index=False)

        log.info(f"alignment stats completed and written to {outdir}/alignment_statistics.xlsx")

    except Exception:

        log.error("Unable to calculate read alignment")
    
     
    log.info("Feature Count from aligned reads started")

    if options.quantification == 'featureCounts':

        fjob = quant.featureCount(gff=options.feature_file, bamDict=outalign, outDir=outdir)

        log.info("feature counts completed and written in %s/Counts.txt",outdir)
    
    elif options.quantification == 'Htseq':

        fjob = quant.htseqCount(gff=options.feature_file, bamDict=outalign, outDir=outdir)


        log.info("feature counts completed and written in %s/Counts.txt",outdir)
    
        # Multi mapped 
    
    if options.mmgg:
        log.info("Counting multimapped read grouops")
       
        mmg_count = mmg.countMMG(sampleDict=samples,bamDict=outalign,gff=options.feature_file+".bed")
        mmg_count.to_excel(os.path.join(outdir,"mmg_count.xlsx", index=False))
        log.info("counting of multimapped read group finished")

    log.info("Differential Expression started with %s",options.detool)

    targets = input_data['targets']
    if options.detool == 'DESeq2':

        count=pd.read_csv(os.path.join(outdir,"Counts_final.txt"), sep="\t")
        result = de.runDESeq2(countDF=count,targetFile=targets,design='sample', combination=combination, subset=False)
        result.to_excel(os.path.join(outdir,"Raw_DEGs_all.xlsx"), index=False)
        
    elif options.detool == 'edgeR':

        count=pd.read_csv(os.path.join(outdir,"Counts_final.txt"), sep="\t")
        result = de.run_edgeR(countDF=count,targetFile=targets, combination=combination, subset=False)
        result.to_excel(os.path.join(outdir,"Raw_DEGs_all.xlsx"), index=False)
        
    log.info("Differential expression analysis completed")
    log.info(f"Filtering differential expressed genes based on logFC {options.fold} and FDR {options.fdr}")
    filtered_DEG= de.degFilter(degDF=result, CompareList=combination,FDR=options.fdr, FOLD=options.fold)
    log.info("filtering DEGs completed ")
    log.info("writting filter DEGs combination wise to excel sheets")
    wd= pd.ExcelWriter(os.path.join(outdir,"filtered_DEGs.xlsx"))
    for key, value in filtered_DEG['filtered'].items():
        value.to_excel(wd,sheet_name=key)
        wd.save()
    wd.close()
    log.info("ploting DEG count figure")
    filtered_DEG['plot'].savefig(os.path.join(outdir,"DEG_figure.png"),dpi=300, bbox_inches='tight')
    # filtered_DEG['plot'].close()
    log.info("Writting DEGs summary to excel file")
    filtered_DEG['summary'].to_excel(os.path.join(outdir,"DEG_count_summary.xlsx"))
    log.info("Creating heatmap of top 50 DEGs")
    if options.heatmap:
        heatmap, ax = pp.plotHeatmap(result,combination,num=50, type=options.heatmaptype)

        heatmap.savefig(os.path.join(outdir,"Top50_gene.png"))

    pu.getGenes(os.path.join(outdir,"filtered_DEGs.xlsx"),combinations=combination, outDir=outdir)

    if options.geneontology:
        outgo = os.path.join(outdir,"Gene_Ontology")
        pu.make_directory(outgo)
        
        for c in combination:
            file = f"{outdir}/diff_genes/{c}.txt"
            ontology_results = go.enrichGO(file =file,species=options.gospecies, type=options.gotype), 
            ontology_results.to_csv(os.path.join(outgo, f"{c}_gene_ontology.txt"), sep="\t", index=False)
    
    if options.keggpathway:
        outkegg = os.path.join(outdir,"KEGG_pathway")
        pu.make_directory(outkegg)
        for c in combination:
            file = f"{outdir}/diff_genes/{c}.txt"
            kegg_results = go.enrichGO(options.keggspecies, file)
            kegg_results.to_csv(os.path.join(outkegg, f"{c}_kegg.txt"), sep="\t", index=False)
    if options.volcanoplot:
        outvolcano = os.path.join(outdir,"Volcano_Plots")
        pu.make_directory(outvolcano)
        for c in combination:
            x,y =pp.plotVolcano(result,c,1)
            x.savefig(outvolcano+"/"+c+"_volcano.png")
        
    if options.maplot:
        outma = os.path.join(outdir,"MA_Plots")
        pu.make_directory(outma)
        for m in combination:
            x,y =pp.plotMA(result,count,m,1)
            x.savefig(outma+"/"+c+"_MA.png")
    
    if options.vennplot:
        degfile = os.path.join(outdir,"filtered_DEGs.xlsx")
        if options.venncombination:
            x = pp.plotVenn(DEGFile=degfile, comparisons=options.venncombination, FOLD=options.fold,outDir=outdir)
            x.savefig(outdir+ "/_Venn.png")
        else:
            if len(combination)<4:
                vnum= len(combination)
            else:
                vnum = len(combination)/4
                
            vlist = np.array_split(combination, vnum)
        
        
            for i in range(0, len(vlist)):
                x = pp.plotVenn(DEGFile=degfile, comparisons=vlist[i], FOLD=options.fold,outDir=outdir)
                x.savefig(outdir+"/Venn_"+i+".png")
        
    endTime = time.ctime()
    log.info("Analysis Complted at %s", endTime)
    log.info("Beer Time!")
        
        
if __name__ == '__main__':
    main()
