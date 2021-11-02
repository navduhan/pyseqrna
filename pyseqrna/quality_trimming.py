#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Title: This modules contains read quality check functions for pySeqRNA
Created : 
@author : Naveen Duhan
'''

import os
import shutil
import sys
import subprocess
import pkg_resources
from pyseqrna.pyseqrna_utils import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu

# Intialize the logger

log = PyseqrnaLogger(mode='a', log='qt')


def flexbarRun(sampleDict,  configFile=None, slurm=False, mem=10, cpu=8, task=1, paired=False, outDir="pySeqRNA_results", dep=''):
    """
    This function is to perform adapter and quality based trimming of reads using flexbar trimming tool

       https://github.com/seqan/flexbar

    Args:
        Samples (Dictionary, required): [Dictionary containing samples]. 
        inputPath (str, required): [input files path]. Defaults to None.
        configFile ([type], required): [Tool parameters file]. Defaults flexbar.ini in settings.
        slurm (bool, optional): [Check if Slurm job scheduling available]. Defaults to False.
        mem (int, required with slurm): [Memory to be used]. Defaults to 10.
        cpu (int, required with slurm): [Number of CPU]. Defaults to 8.
        task (int, required with slurm): [Number of task]. Defaults to 1.
        pairedEND (bool, optional): [True if paired end samples ]. Defaults to False.
        outDir (str, optional): [output directory for results]. Defaults to "pySeqRNA_results".
        dep (str, optional with slurm): [if depends on previsous job]. Defaults to ''.
    """

    if configFile != None:

        try:

            config = pu.parse_config_file(configFile)

        except Exception:

            log.error("Please provide a valid config file")

    else:
        
        log.info("Using default config file flexbar.ini")

        stream = pkg_resources.resource_stream('pyseqrna', "param/flexbar.ini")

        config = pu.parse_config_file(stream.name)

        

    flexbar_config = config[list(config.keys())[0]]

    output = os.path.join(outDir, pu.get_basename(
        str(flexbar_config[0]).split(" ")[1]))

    if os.path.exists(outDir):

        output1 = os.path.join(outDir, output)

        output = pu.make_directory(output1)

    else:
        parent, base = os.path.split(output)
        output = pu.make_directory(base)

    const = ['threads', '-n']  # require when change number of CPUs

    if not slurm:

        flexbar_config = pu.replace_cpu(flexbar_config, const)

    args = ' '.join(flexbar_config[1:])

    outflex = {}  # initialize a dict to carry outfile names for next function

    job_id = []

    for key, value in sampleDict.items():  # Iterate thorough total number of samples

        if paired:
            try:
                input1 = value[2]
                input2 = value[3]

                outflex[key] = [value[0], value[1], os.path.join(
                    output, value[0]+"_1.fastq.gz"), os.path.join(output, value[0]+"_2.fastq.gz")]

            except Exception:

                log.error(
                    "please provide a paired END sample file or input Path is wrong")

            inputPair = pu.get_basename(input1)+" and " + pu.get_basename(input2)

            out = os.path.join(output, value[0])

        else:
            try:
                inputFile = value[2]

            except Exception:

                log.error("please provide a valid input data path")

            outflex[key] = [value[0], value[1],
                            os.path.join(output, value[0]+".fastq.gz")]

            inputPair = pu.getFile(value[2])

            out = os.path.join(output, value[0])

        execPATH = shutil.which("flexbar")

        if execPATH is None:  # get absolute path of flexbar
            
            log.error("flexbar command not found in path")

            sys.exit()

        else:
            
            if paired:

                flexbarCmd = f"{execPATH} -r {input1} -p {input2} -t {out} {args}"

            else:

                flexbarCmd = f"{execPATH} -r {inputFile} -t {out} {args}"

            if slurm == True:  # check if slurm job scheduling is enabled or not

                try:
                    job = pu.clusterRun(job_name='flexbar', sout=os.path.join(output, "flexbar.out"), serror=os.path.join(
                        output, "flexbar.err"), command=flexbarCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)

                    job_id.append(job)

                    log.info(
                        "Job submitted on slurm successfully for {} with {}".format(inputPair, job))

                except Exception:

                    log.error("Slurm job sumission failed")

            else:

                try:
                    with open(os.path.join(output, "flexbar.out"), 'w+') as fout:
                        with open(os.path.join(output, "flexbar.err"), 'w+') as ferr:

                            job = subprocess.call(
                                flexbarCmd, shell=True, stdout=fout, stderr=ferr)

                            job_id.append(" ")

                            log.info(
                                "Job successfully completed for {} with status {}".format(inputPair, job))

                except Exception:

                    log.info("Job submition failed for {} ".format(
                        pu.get_basename(out)))

    return outflex, job_id


def trimmomaticRun(sampleDict=None, configFile=None, slurm=False, mem=10, cpu=8, task=1, paired=False, outDir="pySeqRNA_results", dep=''):
    """
    This function is to perform adapter and quality based trimming of reads using trmmomatic trimming tool


    Args:
        Samples (Dictionary, required): [Dictionary containing samples]. Defaults to None.
        inputPath (str, required): [input files path]. Defaults to None.
        configFile ([type], required): [Tool parameters file]. Defaults flexbar.ini in settings.
        slurm (bool, optional): [Check if Slurm job scheduling available]. Defaults to False.
        mem (int, required with slurm): [Memory to be used]. Defaults to 10.
        cpu (int, required with slurm): [Number of CPU]. Defaults to 8.
        task (int, required with slurm): [Number of task]. Defaults to 1.
        pairedEND (bool, optional): [True if paired end samples ]. Defaults to False.
        outDir (str, optional): [output directory for results]. Defaults to "pySeqRNA_results".
        dep (str, optional with slurm): [if depends on previsous job]. Defaults to ''.
    """

    if configFile != None:

        try:

            config = pu.parse_config_file(configFile)

        except Exception:

            log.error("Please provide a valid config file")

    else:
        if paired:
            log.info("Using default config file trimmomaticPE.ini")
            stream = pkg_resources.resource_stream(
                'pyseqrna', "param/trimmomaticPE.ini")
            config = pu.parse_config_file(stream.name)

        else:
            log.info("Using default config file trimmomaticSE.ini")
            stream = pkg_resources.resource_stream(
                'pyseqrna', "param/trimmomaticSE.ini")
            config = pu.parse_config_file(stream.name)
           

    trimmomatic_config = config[list(config.keys())[0]]

    output = os.path.join(outDir, pu.getFile(str(trimmomatic_config[0])))

    if os.path.exists(outDir):

        output1 = os.path.join(outDir, output)

        output = pu.make_directory(output1)

    else:

        output = pu.make_directory(output)

    const = ['threads', '-n', 'cores']  # require when change number of CPUs

    if not slurm:

        trimmomatic_config = pu.replace_cpu(trimmomatic_config, const)

    args = ' '.join(trimmomatic_config[1:])

    outtrimmomatic = {}  # initialize a dict to carry outfile names for next function

    job_id = []

    for key, value in sampleDict.items():  # Iterate thorough total number of samples

        summ = os.path.join(output, value[0]+"_summary.txt")

        if paired:
            try:
                inputFile1 = value[2]
                inputFile2 = value[3]

                outtrimmomatic[key] = [value[0], value[1], os.path.join(
                    output, value[0]+"_1.fastq.gz"), os.path.join(output, value[0]+"_2.fastq.gz")]

            except Exception:

                log.error(
                    "Please provide a paired END sample file or input Path is wrong")

            out1 = os.path.join(output, value[0]+"_1.fastq.gz")
            out2 = os.path.join(output, value[0]+"_2.fastq.gz")
            inputPair = pu.get_basename(
                inputFile1)+" and " + pu.get_basename(inputFile2)

        else:
            try:
                inputFile = value[2]

            except Exception:

                log.error("please provide a valid input data path")

            outtrimmomatic[key] = [value[0], value[1],
                                   os.path.join(output, value[0]+".fastq.gz")]

            out = os.path.join(output, value[0]+".fastq.gz")

            inputPair = pu.get_basename(inputFile1)

    

        execPATH = shutil.which('trimmomatic') # get absolute path of trimmomatic

        if execPATH is None:

            log.error("trimmomatic command not found in path")

            sys.exit()
        else:

            if paired:

                trimmomaticCmd = f"{execPATH} PE {inputFile1} {inputFile2} {out1} {out2} {args} -summary {summ}"

            else:

                trimmomaticCmd = f"{execPATH} SE {inputFile} {out} {args} -summary {summ}"

            # print(trimmomaticCmd)
            if slurm == True:  # check if slurm job scheduling is enabled or not

                try:
                    job = pu.clusterRun(job_name='trimmomatic', sout=os.path.join(output, "trimmomatic.out"), serror=os.path.join(
                        output, "trimmomatic.err"), command=trimmomaticCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)
                    job_id.append(job)
                    log.info("Job successfully submited for {} with {}".format(
                        inputPair, job))

                except Exception:

                    log.error("Slurm job sumission failed")

            else:

                try:
                    with open(os.path.join(output, "trimmomatic.out"), 'w+') as fout:
                        with open(os.path.join(output, "trimmomatic.err"), 'w+') as ferr:
                            job = subprocess.call(
                                trimmomaticCmd, shell=True, stdout=fout, stderr=ferr)
                            job_id.append(job)
                            log.info(
                                "Job successfully completed for {} ".format(inputPair))

                except Exception:

                    log.error("Job sumission failed")
            try:
                shutil.move("LEADING:3", output+"/LEADING:3")
                shutil.move("TRAILING:3", output+"/TRAILING:3")
            except Exception:
                pass

    return outtrimmomatic, job_id


def trim_galoreRun(sampleDict=None,  configFile=None, slurm=False, mem=10, cpu=8, task=1, paired=False, outDir="pySeqRNA_results", dep=''):
    """
    This function is to perform adapter and quality based trimming of reads using trmmomatic trimming tool


    Args:
        Samples (Dictionary, required): [Dictionary containing samples]. Defaults to None.
        inputPath (str, required): [input files path]. Defaults to None.
        configFile ([type], required): [Tool parameters file]. Defaults flexbar.ini in settings.
        slurm (bool, optional): [Check if Slurm job scheduling available]. Defaults to False.
        mem (int, required with slurm): [Memory to be used]. Defaults to 10.
        cpu (int, required with slurm): [Number of CPU]. Defaults to 8.
        task (int, required with slurm): [Number of task]. Defaults to 1.
        pairedEND (bool, optional): [True if paired end samples ]. Defaults to False.
        outDir (str, optional): [output directory for results]. Defaults to "pySeqRNA_results".
        dep (str, optional with slurm): [if depends on previsous job]. Defaults to ''.
    """

    if configFile != None:

        try:

            config = pu.parse_config_file(configFile)

        except Exception:

            log.error("Please provide a valid config file")

    else:

        stream = pkg_resources.resource_stream(
            'pyseqrna', "param/trim_galore.ini")
        config = pu.parse_config_file(stream.name)
        log.info("Using default config file trim_galore.ini")

    trim_galore_config = config[list(config.keys())[0]]

    output = os.path.join(outDir, pu.getFile(
        str(trim_galore_config[0].split(" ")[1])))

    if os.path.exists(outDir):

        output1 = os.path.join(outDir, output)

        output = pu.make_directory(output1)

    else:

        output = pu.make_directory(output)

    # require when change number of CPUs
    const = ['threads', '-n', 'cores', '-j']

    if not slurm:

        trim_galore_config = pu.replace_cpu(trim_galore_config, const)

    args = ' '.join(trim_galore_config[1:])

    outtrim_galore = {}  # initialize a dict to carry outfile names for next function

    job_id = []

    for key, value in sampleDict.items():  # Iterate thorough total number of samples

        if paired:

            try:
                input1 = value[2]
                input2 = value[3]

                outtrim_galore[key] = [value[0], value[1], os.path.join(
                    output, value[0]+"_val_1.fq.gz"), os.path.join(output, value[0]+"_val_2.fq.gz")]

            except Exception:

                log.error(
                    "please provide a paired END sample file or input Path is wrong")

            out = value[0]
            inputPair = pu.get_basename(
                input1)+" and " + pu.get_basename(input2)

        else:
            try:
                inputFile = value[2]

            except Exception:

                log.error("please provide a valid input data path")

            outtrim_galore[key] = [value[0], value[1],
                                   os.path.join(output, value[0]+"_trimmed.fq.gz")]

            out = value[0]

            inputPair = pu.get_basename(input1)

        

        execPATH = shutil.which('trim_galore')  # get absolute path of trim_galore

        if execPATH is None:

            log.error("trim_galore command not found in path")
            sys.exit()
        else:
            if paired:

                trim_galoreCmd = f"{execPATH} {args} {input1} {input2} --paired --basename {out} -o {output}"

            else:

                trim_galoreCmd = f"{execPATH} {args} {inputFile} --basename {out} -o {output}"

            # print(trim_galoreCmd)

            if slurm == True:  # check if slurm job scheduling is enabled or not

                try:
                    job = pu.clusterRun(job_name='trim_galore', sout=os.path.join(output, "trim_galore.out"), serror=os.path.join(
                        output, "trim_galore.err"), command=trim_galoreCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)
                    job_id.append(job)
                    log.info("Job successfully submited for {} with {}".format(
                        inputPair, job_id))

                except Exception:

                    log.error("Slurm job sumission failed")

            else:

                try:
                    with open(os.path.join(output, "trim_galore.out"), 'w+') as fout:
                        with open(os.path.join(output, "trim_galore.err"), 'w+') as ferr:
                            job = subprocess.call(
                                trim_galoreCmd, shell=True, stdout=fout, stderr=ferr)
                            job_id.append(" ")
                            log.info(
                                "Job successfully completed for {} ".format(inputPair))

                except Exception:

                    log.error("Job sumission failed")

    return outtrim_galore, job_id
