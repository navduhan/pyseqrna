#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
:Title: This module is the main funtion to exceute pySeqRNA

:Created: September 15, 2021

:author: naveen duhan

:version: 0.2
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
from pyseqrna.coexpression import run_Clusts
from pyseqrna import report
from pyseqrna.__organism__ import organism
from tabulate import tabulate
import shutil
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
    
    options = arg_parser.parse_args()  
        
    if options.supported_organism:
        table = [[v, k] for k, v in organism.items()]
            
        print(tabulate(table, headers=["Species Name", "Species code"], tablefmt="grid"))
        
        sys.exit(1)
            
    if options.species is None:
        log.error("Please provide a species and organism type")
        sys.exit()
    if options.species not in organism.keys():
        
        log.info("The provided species does not have functional annotations. Analysis will run up to functional annotation")
        
        
    if options.threads =='80% of available CPU' :
    
        options.threads = pu.get_cpu()
        
    startTime= time.ctime()

    
    # Create directory for results
    if options.resume == 'all':
        dryrun = False
        log.info("Analysis started at %s", startTime)
        outdir = pu.make_directory(options.outdir, dryrun=dryrun)
    else: 
        dryrun = True
        log.info(f"Analysis resume form {options.resume} at {startTime}")
        outdir = pu.make_directory(options.outdir, dryrun=dryrun)
    # with open (os.path.join(outdir,"pyseqrnaa.dill"), 'wb') as dill_save:
        # Read input samples from file
    input_data = pu.read_input_file(options.input_file,options.samples_path, paired = options.paired)

    samples = input_data['samples']
    targets = input_data['targets']
    
    if options.combination != 'all':

        combination = options.combination
    else:

        combination = input_data['combinations']

    if options.resume == 'trimming' or options.resume == 'alignment' or  options.resume == 'differential' or options.resume == 'functional':
        dryrun = True
        qualitydir = pu.make_directory(os.path.join(outdir, "1_Quality"), dryrun=dryrun)
    
        fastqcout = qc.fastqcRun(sampleDict=samples, configFile=options.param, outDir=qualitydir, slurm=options.slurm, mem=options.memory, cpu=options.threads, pairedEND=options.paired, dryrun=dryrun)
    else: 
        log.info("Starting with read quality check")

        dryrun = False
    
        qualitydir = pu.make_directory(os.path.join(outdir, "1_Quality"), dryrun=dryrun)
        
        jobid, fastqcout = qc.fastqcRun(sampleDict=samples, configFile=options.param, outDir=qualitydir, slurm=options.slurm, mem=options.memory, cpu=options.threads, pairedEND=options.paired, dryrun=dryrun)

        if options.slurm:

            for job in jobid:
                wait(lambda: pu.check_status(job), waiting_for="quality to finish")
                log.info(f"Quality check completed for job {job}")

            log.info("Read quality check completed succesfully")
        else:
        
            log.info("Read quality check completed succesfully")

    
    if options.resume == 'alignment' or  options.resume == 'differential' or options.resume == 'functional':
        dryrun = True

        qualitydir = pu.make_directory(os.path.join(outdir, "1_Quality"), dryrun=dryrun)

        if options.trimming == 'flexbar':

            outtrim, jobidt = qt.flexbarRun(sampleDict=samples, configFile=options.param, paired=options.paired, slurm=options.slurm, mem=options.memory,cpu=options.threads, outDir=qualitydir, dryrun=dryrun)
        
        elif options.trimming == 'trimmomatic':

            outtrim, jobidt = qt.trimmomaticRun(sampleDict=samples, configFile=options.param, slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=qualitydir,paired=options.paired, dryrun=dryrun)
            
        elif options.trimming == 'trim_galore':

            outtrim, jobidt = qt.trim_galoreRun(sampleDict=samples, configFile=options.param, slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=qualitydir,paired=options.paired, dryrun=dryrun)
    else: 

        dryrun = False

         # Trimming
        log.info(f"Read trimming started with {options.trimming}")

        if not os.path.isdir(os.path.join(outdir, "1_Quality")):

                qualitydir = pu.make_directory(os.path.join(outdir, "1_Quality"), dryrun=dryrun)

        if options.trimming == 'flexbar':

            outtrim, jobidt = qt.flexbarRun(sampleDict=samples, configFile=options.param, paired=options.paired, slurm=options.slurm, mem=options.memory,cpu=options.threads, outDir=qualitydir, dryrun=dryrun)
        
        elif options.trimming == 'trimmomatic':

            outtrim, jobidt = qt.trimmomaticRun(sampleDict=samples, configFile=options.param, slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=qualitydir,paired=options.paired, dryrun=dryrun)
            
        elif options.trimming == 'trim_galore':

            outtrim, jobidt = qt.trim_galoreRun(sampleDict=samples, configFile=options.param, slurm=options.slurm,mem=options.memory,cpu=options.threads, outDir=qualitydir,paired=options.paired, dryrun=dryrun)
    
        if options.slurm:
            for job in jobidt:
                wait(lambda: pu.check_status(job), waiting_for="trimming to finish")
                log.info(f"Trimming completed for job {job}")

            log.info("Read trimming completed succesfully")
        else:
        
            log.info("Read trimming completed succesfully")

    # Check if read quality check is activated after trimming
    if options.resume == 'alignment' or  options.resume == 'differential' or options.resume == 'functional':
        dryrun = True
        if options.fastqc2:
        
            fastqcout = qc.fastqcRun(sampleDict=outtrim, configFile=options.param, afterTrim=True, outDir=qualitydir, slurm=options.slurm, mem=options.memory, cpu=options.threads, pairedEND=options.paired, dryrun=dryrun)

    else:
        dryrun = False

        if options.fastqc2:
            
            jobid, fastqcout = qc.fastqcRun(sampleDict=outtrim, configFile=options.param, afterTrim=True, outDir=qualitydir, slurm=options.slurm, mem=options.memory, cpu=options.threads, pairedEND=options.paired, dryrun=dryrun)

            if options.slurm:
                for job in jobid:
                    wait(lambda: pu.check_status(job), waiting_for="quality to finish")
                    log.info(f"Quality check for trimmed read completed for job {job}")
        
                log.info("Read quality check for trimmed reads completed succesfully")
            else:
            
                log.info("Read quality check for trimmed reads completed succesfully")

    
    # Removal of ribosomal RNA
    if options.resume == 'alignment' or  options.resume == 'differential' or options.resume == 'functional':
        dryrun = True
        if options.ribosomal:
            outtrim = ribo.sortmernaRun(outtrim, qualitydir, rnaDatabases=options.rnadb, pairedEND= options.paired,cpu=options.threads, slurm=options.slurm, dryrun=dryrun)
    else:
        dryrun = False

        if options.ribosomal:
        

            log.info("Removing ribosomal RNA ")
            
            outtrim, jobs = ribo.sortmernaRun(outtrim, qualitydir, rnaDatabases=options.rnadb, pairedEND= options.paired,cpu=options.threads, slurm=options.slurm, dryrun=dryrun)

            if options.slurm:

                wait(lambda: pu.check_status(jobs), waiting_for="ribosomal rna removal to finish")

                log.info("Ribsomal RNA removal completed succesfully")
            else:
            
                log.info("Ribsomal RNA removal completed succesfully")

    # Alignment Section
    if options.resume == 'differential' or options.resume == 'functional':
        
        dryrun = True

        aligndir = pu.make_directory(os.path.join(outdir, "2_Alignment"), dryrun=dryrun)

        if options.aligner == 'STAR':
            
            aligner= al.STAR_Aligner(genome=options.reference_genome, configFile=options.param,  slurm=options.slurm, memory=options.memory,  outDir=aligndir, dryrun=dryrun)
            jobida = aligner.build_index(mem=options.memory,cpu=options.threads)
            outalign= aligner.run_Alignment(outtrim, pairedEND=options.paired, mem=options.memory, cpu=options.threads)
        
        elif options.aligner == 'hisat2':

            aligner = al.hisat2_Aligner(genome=options.reference_genome, configFile=options.param,  slurm=options.slurm,  outDir=aligndir, dryrun=dryrun)
            jobida= aligner.build_index(mem=options.memory,cpu=options.threads)
            outalign = aligner.run_Alignment(outtrim, pairedEND=options.paired, mem=options.memory, cpu=options.threads)


        quantdir = pu.make_directory(os.path.join(outdir, "3_Quantification"), dryrun=dryrun)

    else: 

        dryrun = False
    
        log.info(f"Starting Alignment  with {options.aligner}")

        aligndir = pu.make_directory(os.path.join(outdir, "2_Alignment"), dryrun=dryrun)

        if options.aligner == 'STAR':

            aligner= al.STAR_Aligner(genome=options.reference_genome, configFile=options.param,  slurm=options.slurm, memory= options.memory, outDir=aligndir, dryrun=dryrun)

            log.info("Genome indexing started")

            jobida = aligner.build_index(mem=options.memory,cpu=options.threads)
    
            log.info("Alignment started")

            outalign, jobalign = aligner.run_Alignment(outtrim, pairedEND=options.paired, mem=options.memory, cpu=options.threads)

        
        elif options.aligner == 'hisat2':

            aligner = al.hisat2_Aligner(genome=options.reference_genome, configFile=options.param,  slurm=options.slurm,  outDir=aligndir, dryrun=dryrun)

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
            align_stat = ps.align_stats(samples,outtrim, outalign,pairedEND=options.paired, cpu=options.threads)

            align_stat.to_excel(aligndir+"/alignment_statistics.xlsx", index=False)

            log.info(f"alignment stats completed and written to {aligndir}/alignment_statistics.xlsx")

        except Exception:

            log.error("Unable to calculate read alignment")
            
        if not options.feature_file.endswith(('gff3', 'gff', 'gtf')):
            
            log.info(" A valid feature file was not provided.")
            
            comb_report = ','.join(combination)
            
            report.generate_report(outdir, comb_report,  options.reference_genome, options.feature_file, options.input_file, options.fold, options.fdr)
                    
            endTime = time.ctime()
            
            log.info("Analysis Complted at %s", endTime)
            
            log.info("Beer Time!")
            
            sys.exit(1)
        
        log.info("Feature Count from aligned reads started")

        quantdir = pu.make_directory(os.path.join(outdir, "3_Quantification"), dryrun=dryrun)

        if options.quantification == 'featureCounts':

            fjob = quant.featureCount(gff=options.feature_file, configFile=options.param, bamDict=outalign, outDir=quantdir, pairedEND=options.paired)

            log.info("feature counts completed and written in %s/Raw_Counts.xlsx",quantdir)
        
        elif options.quantification == 'Htseq':

            fjob = quant.htseqCount(gff=options.feature_file, configFile=options.param, bamDict=outalign, outDir=quantdir)

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
            
        if options.normalizecount == 'FPKM':
            
            r = norm.FPKM()

            r[0].to_excel(os.path.join(quantdir,'FPKM_counts.xlsx'), index=False)

            r[1].savefig(os.path.join(quantdir,'FPKM_vs_Raw_counts.png'), bbox_inches='tight')

        if options.normalizecount == 'CPM':
            
            r = norm.CPM()

            r[0].to_excel(os.path.join(quantdir,'CPM_counts.xlsx'), index=False)

            r[1].savefig(os.path.join(quantdir,'CPM_vs_Raw_counts.png'), bbox_inches='tight')
        
        if options.normalizecount == 'TPM':
            
            r = norm.TPM()

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

    

    if options.resume == 'functional':

        dryrun = True

        diffdir = pu.make_directory(os.path.join(outdir, "4_Differential_Expression"), dryrun= dryrun)
    
    else:

        dryrun = False

        if os.path.exists(os.path.join(outdir, "4_Differential_Expression")):

            shutil.rmtree(os.path.join(outdir, "4_Differential_Expression"))
        
        else:

            diffdir = pu.make_directory(os.path.join(outdir, "4_Differential_Expression"), dryrun= dryrun)


        if options.noreplicate:

            count = pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

            result = de.run_edgeR(countDF=count,targetFile=targets, combination=combination, subset=False, replicate=options.noreplicate)

            ge = de.Gene_Description(species=options.species,combinations=combination, type=options.speciestype, degFile=result, filtered=False)

            results = ge.add_names()

            results.to_excel(os.path.join(diffdir,"All_gene_expression.xlsx"), index=False)
        
        else:

            if options.detool == 'DESeq2':

                count=pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

                result = de.runDESeq2(countDF=count,targetFile=targets,design='sample', combination=combination, subset=False)
                try:
                    ge = de.Gene_Description(species=options.species, type=options.speciestype, degFile=result, filtered=False)

                    results = ge.add_names()

                    results.to_excel(os.path.join(diffdir,"All_gene_expression.xlsx"), index=False)
                    deextraColumns=True
                except Exception:
                    results = result
                    results.to_excel(os.path.join(diffdir,"All_gene_expression.xlsx"), index=False)
                    deextraColumns=False

            elif options.detool == 'edgeR':

                count = pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

                result = de.run_edgeR(countDF=count,targetFile=targets, combination=combination, subset=False, replicate=options.noreplicate)

                ge = de.Gene_Description(species=options.species,combinations=combination, type=options.speciestype, degFile=result, filtered=False)

                results = ge.add_names()

                results.to_excel(os.path.join(diffdir,"All_gene_expression.xlsx"), index=False)
        
        if options.mmgg:

            if options.detool == 'DESeq2':

                countm = pd.read_excel(os.path.join(quantdir,"Raw_MMGcounts.xlsx"))

                resultm = de.runDESeq2(countDF=countm,targetFile=targets,design='sample', combination=combination, subset=False, mmg=True)

                ge = de.Gene_Description(species=options.species, type=options.speciestype, degFile=result, filtered=False)

                resultm.to_excel(os.path.join(diffdir,"All_MMG_expression.xlsx"), index=False)

            elif options.detool == 'edgeR':

                countm = pd.read_excel(os.path.join(quantdir,"Raw_MMGcounts.xlsx"))

                resultm = de.run_edgeR(countDF=countm,targetFile=targets, combination=combination, subset=False, replicate=options.noreplicate, mmg=True)

                ge = de.Gene_Description(species=options.species,combinations=combination, type=options.speciestype, degFile=result, filtered=False)

                resultm.to_excel(os.path.join(diffdir,"All_MMG_expression.xlsx"), index=False)

            
        log.info("Differential expression analysis completed")

        log.info(f"Filtering differential expressed genes based on logFC {options.fold} and FDR {options.fdr}")

        filtered_DEG= de.degFilter(degDF=results, CompareList=combination,FDR=options.fdr, FOLD=options.fold, extraColumns=deextraColumns)

        if options.mmgg:

            filtered_MMG = de.degFilter(degDF=resultm, CompareList=combination, FDR=options.fdr, FOLD=options.fold, mmg=True, extraColumns=False)

        log.info("filtering DEGs completed ")

        log.info("writting filter DEGs combination wise to excel sheets")
        # write up and down filtered genes together
        wa = pd.ExcelWriter(os.path.join(diffdir,"Filtered_DEGs.xlsx"))

        for key, value in filtered_DEG['filtered'].items():
            value.to_excel(wa,sheet_name=key, index=False)
            
        wa.close()

        # write up filtered genes together
        wu = pd.ExcelWriter(os.path.join(diffdir,"Filtered_upDEGs.xlsx"))

        for key, value in filtered_DEG['filteredup'].items():
            value.to_excel(wu,sheet_name=key, index=False)
            
        wu.close()
        # write down filtered genes together
        wd = pd.ExcelWriter(os.path.join(diffdir,"Filtered_downDEGs.xlsx"))

        for key, value in filtered_DEG['filtereddown'].items():
            value.to_excel(wd,sheet_name=key, index=False)
            
        wd.close()

        if options.mmgg:
            wa = pd.ExcelWriter(os.path.join(diffdir,"Filtered_MMGs.xlsx"))

            for key, value in filtered_MMG['filtered'].items():
                value.to_excel(wa,sheet_name=key, index=False)
                
            wa.close()

            # write up filtered genes together
            wu = pd.ExcelWriter(os.path.join(diffdir,"Filtered_upMMGs.xlsx"))

            for key, value in filtered_MMG['filteredup'].items():
                value.to_excel(wu,sheet_name=key, index=False)
                
            wu.close()
            # write down filtered genes together
            wd = pd.ExcelWriter(os.path.join(diffdir,"Filtered_downMMGs.xlsx"))

            for key, value in filtered_MMG['filtereddown'].items():
                value.to_excel(wd,sheet_name=key, index=False)
                
            wd.close()

        log.info("ploting DEG count figure")

        filtered_DEG['plot'].savefig(os.path.join(diffdir,"Filtered_DEG.png"),dpi=300, bbox_inches='tight')

        if options.mmgg:

            log.info("ploting MMG count figure")

            filtered_MMG['plot'].savefig(os.path.join(diffdir,"Filtered_MMG.png"),dpi=300, bbox_inches='tight')

        log.info("Writting DEGs summary to excel file")
        
        filtered_DEG['summary'].to_excel(os.path.join(diffdir,"Filtered_DEGs_summary.xlsx"), index=False)

        if options.mmgg:
            log.info("Writting MMGs summary to excel file")
        
            filtered_MMG['summary'].to_excel(os.path.join(diffdir,"Filtered_MMGs_summary.xlsx"), index=False)

        genesdir = pu.getGenes(os.path.join(diffdir,"Filtered_DEGs.xlsx"),combinations=combination, outDir=diffdir)

        if options.mmgg:

            genesdir = pu.getGenes(os.path.join(diffdir,"Filtered_MMGs.xlsx"),combinations=combination, outDir=diffdir, mmg=True)

        if options.coexpression:

            coexpout = pu.make_directory(os.path.join(diffdir, 'co-expression'))

            run_Clusts(countFile=os.path.join(quantdir,"Raw_Counts.xlsx"), targets=targets, outdir=coexpout)

            


    if  options.resume == 'functional' or options.resume == 'all':

        if os.path.exists(os.path.join(outdir, "5_Visualization")):

            plotdir = os.path.join(outdir, "5_Visualization")

        else:
    
            plotdir = pu.make_directory(os.path.join(outdir, "5_Visualization"))

            
        if options.heatmap:

            log.info("Creating heatmap of top 50 DEGs")
            
            if options.heatmaptype == 'degs':
                
                degs = pd.read_excel(os.path.join(diffdir,"All_gene_expression.xlsx"))

                heatmap, ax = pp.plotHeatmap(degs,combination,num=50, type=options.heatmaptype)

                heatmap.savefig(os.path.join(plotdir,f"Heatmap_top50.png"), bbox_inches='tight')
                
            if options.heatmaptype == 'counts':
                
                count = pd.read_excel(os.path.join(quantdir,"Raw_Counts.xlsx"))

                heatmap, ax = pp.plotHeatmap(count,combination,num=50, type=options.heatmaptype)

                heatmap.savefig(os.path.join(plotdir,f"Heatmap_top50.png"), bbox_inches='tight')
                
        if options.pcaplot:
            
            log.info(f"creating PCA plot based on {options.pcaplottype} values")
            
            norm = Normalization(countFile= os.path.join(quantdir,"Raw_Counts.xlsx"), featureFile=options.feature_file, keyType=options.source)
            
            if options.pcaplottype == 'RPKM':
                
                r = norm.RPKM()

                pcafig = pp.pcaPlot(ncountdf=r[0])
                pcafig.savefig(os.path.join(plotdir,f"PCA_plot.png"), bbox_inches='tight')
                
            if options.pcaplottype == 'FPKM':
                
                r = norm.FPKM()

                pcafig = pp.pcaPlot(ncountdf=r[0])
                pcafig.savefig(os.path.join(plotdir,f"PCA_plot.png"), bbox_inches='tight')

            if options.pcaplottype == 'CPM':
                
                r = norm.CPM()

                pcafig = pp.pcaPlot(ncountdf=r[0])
                pcafig.savefig(os.path.join(plotdir,f"PCA_plot.png"), bbox_inches='tight')
                
            if options.pcaplottype == 'TPM':
                
                r = norm.TPM()
                
                pcafig = pp.pcaPlot(ncountdf=r[0])
                pcafig.savefig(os.path.join(plotdir,f"PCA_plot.png"), bbox_inches='tight')
                
            if options.pcaplottype == 'medianRatiocount':
                
                r = norm.meanRatioCount()

                pcafig = pp.pcaPlot(ncountdf=r)
                pcafig.savefig(os.path.join(plotdir,f"PCA_plot.png"), bbox_inches='tight')
            
        if options.species not in organism.keys():
        
            log.info("The provided species is nbot currently supported for functional annotations.")  
        
            comb_report = ','.join(combination)
            report.generate_report(outdir, comb_report,  options.reference_genome, options.feature_file, options.input_file, options.fold, options.fdr)
                    
            endTime = time.ctime()
            log.info("Analysis Complted at %s", endTime)
            log.info("Beer Time!")

            sys.exit(1)
    
        if os.path.exists(os.path.join(outdir, "6_Functional_Annotation")):

            annodir = os.path.join(outdir, "6_Functional_Annotation")

        else:

            annodir = pu.make_directory(os.path.join(outdir, "6_Functional_Annotation"))

        if options.geneontology:

            if os.path.exists(os.path.join(annodir,"Gene_Ontology")):

                shutil.rmtree(os.path.join(annodir,"Gene_Ontology"))

            else:

                outgo = pu.make_directory(os.path.join(annodir,"Gene_Ontology"))

                gofiles = pu.make_directory(os.path.join(outgo,"GO_Files"))

                goplots = pu.make_directory(os.path.join(outgo,"GO_Plots"))
            
                if options.mmgg:

                    gofilesM = pu.make_directory(os.path.join(outgo,"GO_Files_MMG"))

                    goplotsM = pu.make_directory(os.path.join(outgo,"GO_Plots_MMG"))

            go = GeneOntology(species=options.species, type=options.speciestype, keyType=options.source,  gff=options.feature_file)

            for c in combination:

                file_deg = f"{genesdir}/{c}.txt"

                ontology_results = go.enrichGO(file=file_deg)

                if ontology_results != "No Gene Ontology":
                    try:
                        go_results = ge.add_names_annotation(ontology_results['result'])

                        go_results.to_excel(f"{gofiles}/{c}_gene_ontology.xlsx", index=False)

                        ontology_results['plot'].savefig(f"{goplots}/{c}_go_dotplot.png", bbox_inches='tight')
                        plt.close()
                    except:
                        continue
                else:
                    log.info(f"No ontology found in {c}")

            if options.mmgg:

                for c in combination:

                    file_mmg = f"{genesdir}/{c}_mmg.txt"

                    ontology_results_mmg = go.enrichGO(file=file_mmg)

                    if ontology_results_mmg != "No Gene Ontology":
                        
                        go_results_MMG = ge.add_names_annotation(ontology_results_mmg['result'])

                        go_results_MMG_go = pu.add_MMG(degDF=os.path.join(diffdir,"Filtered_MMGs.xlsx"), anotDF=go_results_MMG, combination=c)

                        go_results_MMG_go.to_excel(f"{gofilesM}/{c}_mmg_gene_ontology.xlsx", index=False)

                        ontology_results_mmg['plot'].savefig(f"{goplotsM}/{c}_mmg_go_dotplot.png", bbox_inches='tight')

                        plt.close()

                    else:
                        
                        log.info(f"No ontology found in {c}")
        
        if options.keggpathway:


            if os.path.exists(os.path.join(annodir,"KEGG_Pathway")):

                shutil.rmtree(os.path.join(annodir,"KEGG_Pathway"))
            
            else:

                outkegg = pu.make_directory(os.path.join(annodir,"KEGG_Pathway"))

                keggfiles = pu.make_directory(os.path.join(outkegg,"KEGG_Files"))

                keggplots = pu.make_directory(os.path.join(outkegg,"KEGG_Plots"))

                if options.mmgg:

                    keggfilesM = pu.make_directory(os.path.join(outkegg,"KEGG_Files_MMG"))

                    keggplotsM = pu.make_directory(os.path.join(outkegg,"KEGG_Plots_MMG"))

            pt = Pathway(species=options.species, keyType=options.source, gff= options.feature_file)

            for c in combination:

                file_deg = f"{genesdir}/{c}.txt"

                kegg_results = pt.enrichKEGG(file_deg)

                if kegg_results != "No Pathway":
                    try:
                        k_result = ge.add_names_annotation(kegg_results['result'])

                        k_result.to_excel(os.path.join(keggfiles, f"{c}_kegg.xlsx"), index=False)

                        kegg_results['plot'].savefig( f"{keggplots}/{c}_kegg_dotplot.png", bbox_inches='tight')

                        plt.close()
                    except:
                        continue

                else:
                    log.info(f"No pathway found in {c}")

            if options.mmgg:

                for c in combination:

                    file_mmgs = f"{genesdir}/{c}_mmg.txt"

                    kegg_results_mmg = pt.enrichKEGG(file=file_mmgs)

                    if kegg_results_mmg != "No Pathway":

                        kegg_results_MMG = ge.add_names_annotation(ontology_results_mmg['result'])

                        kegg_results_MMG_go = pu.add_MMG(degDF=os.path.join(diffdir,"Filtered_MMGs.xlsx"), anotDF=kegg_results_MMG, combination=c)
                        
                        kegg_results_MMG_go.to_excel(f"{keggfilesM}/{c}_mmg_kegg.xlsx", index=False)

                        kegg_results_mmg['plot'].savefig(f"{keggplotsM}/{c}_mmg_kegg_dotplot.png", bbox_inches='tight')
                        
                        plt.close()
                        
                    else:
                        log.info(f"No ontology found in {c}")

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
    report.generate_report(outdir, comb_report,  options.reference_genome, options.feature_file, options.input_file, options.fold, options.fdr)
            
    endTime = time.ctime()
    log.info("Analysis Complted at %s", endTime)
    log.info("Beer Time!")
        
        
if __name__ == '__main__':
    main()
