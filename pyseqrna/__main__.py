#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
:Title: This module is the main funtion to exceute pySeqRNA

:Created: September 15, 2021

:author: naveen duhan

:version: 0.1
"""

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
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
from pyseqrna import clustering as cl
from pyseqrna import multimapped_groups as mmg
from pyseqrna.gene_ontology import GeneOntology
from pyseqrna.pathway import Pathway
from pyseqrna.normalize_counts import Normalization
from pyseqrna import report

import math
import pyseqrna.version
import pandas as pd
import numpy as np
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
    # with open (os.path.join(outdir,"pyseqrnaa.dill"), 'wb') as dill_save:
        # Read input samples from file
    input_data = pu.read_input_file(options.input_file,options.samples_path, paired = options.paired)

    samples = input_data['samples']
    targets = input_data['targets']
    
    if options.combination != 'all':

        combination = options.combination
    else:

        combination = input_data['combinations']

    log.info("Starting with read quality check")

    
    qualitydir = pu.make_directory(os.path.join(outdir, "1_Quality"))
    
    jobid, fastqcout = qc.fastqcRun(sampleDict=samples,outDir=qualitydir, slurm=options.slurm, mem=options.memory, cpu=options.threads, pairedEND=options.paired)

    if options.slurm:
        for job in jobid:
            wait(lambda: pu.check_status(job), waiting_for="quality to finish")
            log.info(f"Quality check completed for job {job}")

        log.info("Read quality check completed succesfully")
    else:
    
        log.info("Read quality check completed succesfully")

    
    # Trimming

    log.info(f"Read trimming started with {options.trimming}")


    if not os.path.isdir(os.path.join(outdir, "1_Quality")):

            qualitydir = pu.make_directory(os.path.join(outdir, "1_Quality"))

    if options.trimming == 'flexbar':

        outtrim, jobidt = qt.flexbarRun(sampleDict=samples,paired=options.paired, slurm=options.slurm, mem=options.memory,cpu=options.threads, outDir=qualitydir)
    
    elif options.trimming == 'trimmomatic':

        outtrim, jobidt = qt.trimmomaticRun(sampleDict=samples,slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=qualitydir,paired=options.paired)
        
    elif options.trimming == 'trim_galore':

        outtrim, jobidt = qt.trim_galoreRun(sampleDict=samples,slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=qualitydir,paired=options.paired)
    
    # dill.dump([input_data, outtrim],dill_save,protocol=dill.HIGHEST_PROTOCOL)
    if options.slurm:
        for job in jobidt:
            wait(lambda: pu.check_status(job), waiting_for="trimming to finish")
            log.info(f"Trimming completed for job {job}")

        log.info("Read trimming completed succesfully")
    else:
    
        log.info("Read trimming completed succesfully")

    # Check if read quality check is activated after trimming
    
    if options.fastqc2:
        
        jobid, fastqcout = qc.fastqcRun(sampleDict=outtrim,afterTrim=True, outDir=qualitydir, slurm=options.slurm, mem=options.memory, cpu=options.threads, pairedEND=options.paired)

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
        
        outtrim, jobs = ribo.sortmernaRun(outtrim, qualitydir, rnaDatabases=options.rnadb, pairedEND= options.paired,cpu=options.threads, slurm=options.slurm)

        if options.slurm:

            wait(lambda: pu.check_status(jobs), waiting_for="ribosomal rna removal to finish")

            log.info("Ribsomal RNA removal completed succesfully")
        else:
        
            log.info("Ribsomal RNA removal completed succesfully")

    # Alignment Section
    
    log.info("Starting Alignment process")

    aligndir = pu.make_directory(os.path.join(outdir, "2_Alignment"))

    if options.aligners == 'STAR':

        aligner= al.STAR_Aligner(genome=options.reference_genome,  slurm=options.slurm,  outDir=aligndir)

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

        aligner = al.hisat2_Aligner(genome=options.reference_genome,  slurm=options.slurm,  outDir=aligndir)

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

        align_stat.to_excel(aligndir+"/alignment_statistics.xlsx", index=False)

        log.info(f"alignment stats completed and written to {aligndir}/alignment_statistics.xlsx")

    except Exception:

        log.error("Unable to calculate read alignment")
    
    log.info("Feature Count from aligned reads started")

    quantdir = pu.make_directory(os.path.join(outdir, "3_Quantification"))

    if options.quantification == 'featureCounts':

        fjob = quant.featureCount(gff=options.feature_file, bamDict=outalign, outDir=quantdir)

        log.info("feature counts completed and written in %s/Raw_Counts.xlsx",quantdir)
    
    elif options.quantification == 'Htseq':

        fjob = quant.htseqCount(gff=options.feature_file, bamDict=outalign, outDir=quantdir)

        log.info("feature counts completed and written in %s/Raw_Counts.txt",quantdir)
    
        # Multi mapped 
    
    if options.mmgg:
        log.info("Counting multimapped read groups")
    
        mmg_count = mmg.countMMG(sampleDict=samples,bamDict=outalign,gff=options.feature_file, minCount=options.minreadcount, percentSample=options.percentsample)

        mmg_count.to_excel(os.path.join(quantdir,"Raw_MMGcounts.xlsx" ),index=False)

        log.info("counting of multimapped read group finished")

    log.info("Differential Expression started with %s",options.detool)

    log.info("Converting Raw Counts to normalized counts")


    norm = Normalization(countFile= os.path.join(quantdir,"Raw_Counts.xlsx"), featureFile=options.feature_file, keyType=options.source)

    if options.normalizecount == 'RPKM':
        
        r = norm.RPKM()

        r[0].to_excel(os.path.join(quantdir,'RPKM_counts.xlsx'), index=False)

        r[1].savefig(os.path.join(quantdir,'RPKM_vs_Raw_counts.png'), bbox_inches='tight')

    if options.normalizecount == 'CPM':
        
        r = norm.CPM()

        r[0].to_excel(os.path.join(quantdir,'CPM_counts.xlsx'), index=False)

        r[1].savefig(os.path.join(quantdir,'CPM_vs_Raw_counts.png'), bbox_inches='tight')
    
    if options.normalizecount == 'TPM':
        
        r = norm.CPM()

        r[0].to_excel(os.path.join(quantdir,'TPM_counts.xlsx'), index=False)

        r[1].savefig(os.path.join(quantdir,'TPM_vs_Raw_counts.png'), bbox_inches='tight')

    if options.normalizecount == 'medianRatiocount':
        
        r = norm.meanRatioCount()

        r.to_excel(os.path.join(quantdir,'Median_ratio_counts.xlsx'),index=False)

    log.info("Clustering samples based on similarity")

    if options.cluster:
        count = pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))
        plot = cl.clusterSample(countDF= count)
        plot.savefig(os.path.join(quantdir,'Sample_cluster.png'), bbox_inches='tight')

    targets = input_data['targets']

    diffdir = pu.make_directory(os.path.join(outdir, "4_Differential_Expression"))


    if options.detool == 'DESeq2':

        count=pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

        result = de.runDESeq2(countDF=count,targetFile=targets,design='sample', combination=combination, subset=False)

        ge = de.Gene_Description(species=options.species, type=options.speciestype, degFile=result, filtered=False)

        results = ge.add_names()

        results.to_excel(os.path.join(diffdir,"All_gene_expression.xlsx"), index=False)

        
    elif options.detool == 'edgeR':

        count = pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

        result = de.run_edgeR(countDF=count,targetFile=targets, combination=combination, subset=False, replicate=options.noreplicate)

        ge = de.Gene_Description(species=options.species,combinations=combination, type=options.speciestype, degFile=result, filtered=False)

        results = ge.add_names()

        results.to_excel(os.path.join(diffdir,"All_gene_expression.xlsx"), index=False)
        
    log.info("Differential expression analysis completed")

    log.info(f"Filtering differential expressed genes based on logFC {options.fold} and FDR {options.fdr}")

    filtered_DEG= de.degFilter(degDF=results, CompareList=combination,FDR=options.fdr, FOLD=options.fold, extraColumns=True)

    log.info("filtering DEGs completed ")

    log.info("writting filter DEGs combination wise to excel sheets")
    # write up and down filtered genes together
    wa = pd.ExcelWriter(os.path.join(diffdir,"Filtered_DEGs.xlsx"))

    for key, value in filtered_DEG['filtered'].items():
        value.to_excel(wa,sheet_name=key, index=False)
        wa.save()
    wa.close()

    # write up filtered genes together
    wu = pd.ExcelWriter(os.path.join(diffdir,"Filtered_upDEGs.xlsx"))

    for key, value in filtered_DEG['filteredup'].items():
        value.to_excel(wu,sheet_name=key, index=False)
        wu.save()
    wu.close()
    # write down filtered genes together
    wd = pd.ExcelWriter(os.path.join(diffdir,"Filtered_downDEGs.xlsx"))

    for key, value in filtered_DEG['filtereddown'].items():
        value.to_excel(wd,sheet_name=key, index=False)
        wd.save()
    wd.close()

    log.info("ploting DEG count figure")

    filtered_DEG['plot'].savefig(os.path.join(diffdir,"Filtered_DEG.png"),dpi=300, bbox_inches='tight')

    log.info("Writting DEGs summary to excel file")

    filtered_DEG['summary'].to_excel(os.path.join(diffdir,"Filtered_DEGs_summary.xlsx"), index=False)

    plotdir = pu.make_directory(os.path.join(outdir, "5_Visualization"))

    if options.heatmap:

        log.info("Creating heatmap of top 50 DEGs")

        heatmap, ax = pp.plotHeatmap(result,combination,num=50, type=options.heatmaptype)

        heatmap.savefig(os.path.join(plotdir,f"Heatmap_top50.png"), bbox_inches='tight')

    genesdir = pu.getGenes(os.path.join(diffdir,"Filtered_DEGs.xlsx"),combinations=combination, outDir=diffdir)

    annodir = pu.make_directory(os.path.join(outdir, "6_Functional_Annotation"))

    if options.geneontology:

        outgo = pu.make_directory(os.path.join(annodir,"Gene_Ontology"))

        gofiles = pu.make_directory(os.path.join(outgo,"GO_Files"))

        goplots = pu.make_directory(os.path.join(outgo,"GO_Plots"))

        go = GeneOntology(species=options.species, type=options.speciestype, keyType=options.source, taxid=options.taxid, gff=options.feature_file)

        for c in combination:

            file_deg = f"{genesdir}/{c}.txt"

            ontology_results = go.enrichGO(file=file_deg)

            if ontology_results != "No Gene Ontology":

                go_results = ge.add_names_annotation(ontology_results['result'])

                go_results.to_excel(f"{gofiles}/{c}_gene_ontology.xlsx", index=False)

                ontology_results['plot'].savefig(f"{goplots}/{c}_go_dotplot.png", bbox_inches='tight')
                plt.close()
            else:
                log.info(f"No ontology found in {c}")
    
    if options.keggpathway:

        outkegg = pu.make_directory(os.path.join(annodir,"KEGG_Pathway"))

        keggfiles = pu.make_directory(os.path.join(outkegg,"KEGG_Files"))

        keggplots = pu.make_directory(os.path.join(outkegg,"KEGG_Plots"))

        pt = Pathway(species=options.species, keyType=options.source, gff= options.feature_file)

        for c in combination:

            file_deg = f"{genesdir}/{c}.txt"

            kegg_results = pt.enrichKEGG(file_deg)

            if kegg_results != "No Pathway":

                k_result = ge.add_names_annotation(kegg_results['result'])

                k_result.to_excel(os.path.join(keggfiles, f"{c}_kegg.xlsx"), index=False)

                kegg_results['plot'].savefig( f"{keggplots}/{c}_kegg_dotplot.png", bbox_inches='tight')

                plt.close()

            else:
                log.info(f"No pathway found in {c}")

    if options.volcanoplot:

        outvolcano = os.path.join(plotdir,"Volcano_Plots")

        pu.make_directory(outvolcano)

        for c in combination:

            x = pp.plotVolcano(result,c,FOLD=options.fold)

            if type(x) != str:

                x.savefig(outvolcano+"/"+c+"_volcano.png")
                plt.close()

    if options.maplot:

        count = pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

        outma = pu.make_directory(os.path.join(plotdir,"MA_Plots"))

        for m in combination:

            x = pp.plotMA(degDF= result,countDF= count,comp=m,FOLD=options.fold, FDR=options.fdr)

            if type(x) != str:

                x.savefig(outma+"/"+m+"_MA.png")
                plt.close()

    if options.vennplot:

        outvenn = pu.make_directory(os.path.join(plotdir,"Venn_Plots"))
        
        degfile = os.path.join(diffdir,"Filtered_DEGs.xlsx")

        if options.venncombination != 'random':

            x = pp.plotVenn(DEGFile=degfile, comparisons=options.venncombination, FOLD=options.fold)

            x.savefig(f"{outvenn}/Venn0.png")

            plt.close()

        else:
            if len(combination)<4:

                vnum = len(combination)

            else:

                vnum = math.ceil(len(combination)/4)
                
            vlist = np.array_split(combination, vnum)
        
            for i in range(0, len(vlist)):

                x = pp.plotVenn(DEGFile=degfile, comparisons= vlist[i], FOLD=options.fold, degLabel=None)

                x.savefig(outvenn+"/Venn_"+str(i)+".png")

                plt.close()
    
    comb_report = ','.join(combination)
    report.generate_report(outdir, comb_report, options.input_file, options.fold, options.fdr)
            
    endTime = time.ctime()
    log.info("Analysis Complted at %s", endTime)
    log.info("Beer Time!")
        
        
if __name__ == '__main__':
    main()
