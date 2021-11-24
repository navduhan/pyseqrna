#!/usr/bin/env python
'''
Title: This script contains quality trimming tools for pySeqRNA
Author: Naveen Duhan
Version: 0.1
'''

import pandas as pd 
import pysam
import pyfastx
import pandas as pd
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu
import matplotlib as plt
import multiprocessing

log = PyseqrnaLogger(mode='a', log="stats")


def getNreads(file, rdict, sp):
    """
        Get total number of reads in fastq file
    Args:
        file ([type]): fastq file

    Returns:
        [type]: total number of reads
    """
    result = len(pyfastx.Fastq(file))
    
    rdict[sp] = int(result)
    # result = subprocess.check_output(f"gzcat {file} | echo $((`wc -l`/4))", shell=True).decode('utf-8').rstrip()
    # logger.info(f"{result} input reads in {file}")
    
    return rdict

def getAligned_reads(file, rdict, sp):

    aligned=0

    for read in pysam.AlignmentFile(file,'rb'):
        if read.is_unmapped == False and read.is_secondary==False:
            aligned += 1
    rdict[sp] = int(aligned)
    return rdict


def getUniquely_mapped(file, rdict, sp):

    uniq_mapped = 0
    for read in pysam.AlignmentFile(file,'rb'):
    
        try:
            if read.get_tag('NH')==1 and read.is_secondary==False:
                uniq_mapped += 1
        except:
            pass
    rdict[sp] = int(uniq_mapped)
    return rdict
    

def getMulti_mapped(file,rdict, sp):

    multi_mapped = 0 

    for read in pysam.AlignmentFile(file,'rb'):
        try:
            if read.get_tag('NH') > 1 and read.is_secondary==False:
                multi_mapped += 1
        except:
            pass
    rdict[sp] = int(multi_mapped)
    return rdict
    


def align_stats(sampleDict=None,trimDict=None, bamDict=None,ribodict=None, pairedEND=False):
        
    manager = multiprocessing.Manager()
    Ireads = manager.dict()
    Nreads = manager.dict()
    Rreads = manager.dict()
    Areads = manager.dict()
    Ureads = manager.dict()
    Mreads = manager.dict()
    processes = []
    
    for sp in sampleDict:
        try:
            p=multiprocessing.Process(target= getNreads, args=(sampleDict[sp][2],Ireads, sp))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
       
            if pairedEND:
                Ireads.update((x, y*2) for x, y in Ireads.items())
            
        except Exception:
            log.error(f"Not able to count Input read number in {sp}")
   
    
    for tf in trimDict:
        try:
            p=multiprocessing.Process(target= getNreads, args=(trimDict[tf][2],Nreads, tf))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
       
            if pairedEND:
                Nreads.update((x, y*2) for x, y in Nreads.items())
        except Exception:
            log.error(f"Not able to count Trim read number in {tf}")

    for bf in bamDict:
        p=multiprocessing.Process(target= pysam.sort, args=(bamDict[bf][2], "-o", bamDict[bf][2].split(".bam")[0] + "_sorted.bam"))
        
        processes.append(p)
        p.start()
            
        for process in processes:
            process.join()
        # pysam.sort(bamDict[bf][2], "-o", bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
    for bf in bamDict:   
        file = bamDict[bf][2].split(".bam")[0] + "_sorted.bam"
        print(file)
        p=multiprocessing.Process(target= pysam.index, args=(file))
        
        processes.append(p)
        p.start()
            
        for process in processes:
            process.join()
        # pysam.index(bamDict[bf][2].split(".bam")[0] + "_sorted.bam")
    for bf in bamDict:
        try:
            p=multiprocessing.Process(target= getAligned_reads, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Areads, bf))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
       
        except Exception:
            log.error(f"Not able to count Aligned read number in {bf}")
        try:
            p=multiprocessing.Process(target= getUniquely_mapped, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Ureads, bf))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()
                
        except Exception:
            log.error(f"Not able to count Uniquely mapped read number in {bf}")
        try:
            p=multiprocessing.Process(target= getMulti_mapped, args=(bamDict[bf][2].split(".bam")[0] + "_sorted.bam",Mreads, bf))
        
            processes.append(p)
            p.start()
            
            for process in processes:
                process.join()           
        except Exception:
            log.error(f"Not able to count Multi mapped read number in {bf}")
    try:
        total = []
        for k in Ireads:

            total.append([k,Ireads[k],Nreads[k],round(Nreads[k]/Ireads[k]*100,2),
                    Areads[k],round(Areads[k]/Nreads[k]*100,2),Ureads[k],round(Ureads[k]/Areads[k]*100,2),Mreads[k], round(Mreads[k]/Areads[k]*100,2)])
        if pairedEND:

            totalDF = pd.DataFrame(total,columns=['Sample','Input_reads2x','Cleaned2x', '%_Cleaned2x','Aligned','%_Aligned','Uniquely_mapped','%_Uniquely_mapped','Multi_mapped', '%_Multi_mapped'])
        else:
            totalDF = pd.DataFrame(total,columns=['Sample','Input_reads','Cleaned', '%_Cleaned','Aligned','%_Aligned','Uniquely_mapped','%_Uniquely_mapped','Multi_mapped', '%_Multi_mapped'])
    except Exception:
        log.error(f"Not able to generate align stats")


    return totalDF


