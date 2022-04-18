#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Title: This modules contains utility functions for pySeqRNA
Created : 
@author : Naveen Duhan
'''



import os
import re 
import sys
import math
import fnmatch
import psutil
import logging
import subprocess
import configparser
import pandas as pd


def PyseqrnaLogger(mode, log):

    """[intialize logger in the pySeqRNA modules]

    Args:
        logger (str): [logger name for the module]
        logFile (str): [file name for logging]
    """
    logger = logging.getLogger(log)
    logger.propagate=False

    # set format for logging
    logFormatter =logging.Formatter('[%(asctime)s]  %(module)s :: %(levelname)s : %(message)s',datefmt='%H:%M:%S')
    # logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # write log in a user provided log file
    
    fileHandler = logging.FileHandler("{}".format('pyseqrna.log'), mode= mode)

    fileHandler.setFormatter(logFormatter)

    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()

    consoleHandler.setFormatter(logging.Formatter('[%(asctime)s]  %(module)s :: %(levelname)s : %(message)s',datefmt='%H:%M:%S'))

    consoleHandler.setLevel(logging.DEBUG)

    logger.addHandler(consoleHandler)

    return logger

log = PyseqrnaLogger(mode='a', log='pu')

def read_input_file(infile, inpath, paired = False):

    """This function reads input sample file and convert into a dictionary. 
        It also make all possible combination for DEG analysis.
        Target dataframe for differential analysis.

    Args:
        inputFile (str, required): [input sample file containing the infromation about project].
        inputPath ([type], required): [Path for input fastq files]. 
        pairedEND (bool, required): [Check if reads are paired end]. Defaults to False.

    Returns:
        dict: [Contains samples, combinations and target]
    """
    samples = {} # dictionary to collect sample information from input sample file
    factors = [] # list for collecting the Identifier infromation from input sample file
    combinations = [] # list for genrating combinations for differential expression analysis

    try:
        with open(infile) as file:

            log.info("Reading input samples File ")

            for line in file:

                if not line.startswith("#") and not line.startswith("SampleName"):
                    line = line.strip()
                    lines = re.split('\s+', line.rstrip())

                    if paired:

                        samples[lines[1]] = [lines[1],lines[2],os.path.join(inpath,lines[3]), os.path.join(inpath,lines[4])]

                    else:

                        samples[lines[1]] = [lines[1], lines[2], os.path.join(inpath,lines[3])]

                    if lines[2] not in factors:

                        factors.append(lines[2])
    except Exception:

        log.error("Please provide a valid input file")

        sys.exit()

    finally:

        log.info("Input file %s read succesfully", infile)

    try:

        # create combinations from factors 
        for i in factors:
            for j in factors:

                if i != j:

                    if j+"-"+i not in combinations:

                        combinations.append(i+"-"+j)
    except Exception:

        log.error("Please provide a valid input file")
    
    finally:

        log.info("Combination created succesfully from %s", infile)

        samplename = []
        sample = []

    try:
        for k, s in samples.items():

            samplename.append(k)

            sample.append(s[1])

        targets = pd.DataFrame(samplename,index=[i for i in samplename])

        targets = targets.assign(sample=sample)

    except Exception:

        log.error("Please provide a valid input file")
    
    finally:

        log.info("targets dataframe for differenatial created succesfully from %s", infile)

    return {'samples': samples, "combinations": combinations, "targets": targets}


def parse_config_file(infile):
    """
    This function parse the config file for all the programs used in pySeqRNA

    Args:
        configFile ([type], optional): [description]. 

    retrun:
        rtype: a dictionary

    """
    sections_dict = {}

    config = configparser.ConfigParser()

    try:

        config.read([infile])

        sections = config.sections()

        for section in sections:

            options = config.options(section)

            temp_dict = {}
            voption = []

            for option in options:

                cc = config.get(section, option)

                temp_dict[option] = cc

            for k, value in temp_dict.items():

                if 'NA' not in value:

                    voption.append(value)

            sections_dict[section] = voption

    except Exception:

        log.error("Please provide a valid config file")
    
    finally:

        log.info("Config generated succesfully from %s", infile)

    return sections_dict

def clusterRun(job_name='pyseqRNA',sout=" pyseqrna", serror="pyseqrna", command='command', time=4, mem=10, cpu=8, tasks=1, dep=''):
    """
    This function is for submitting job on cluster with SLURM job Scheduling

    Args:
        job_name (str, optional): [description]. Defaults to 'pyseqRNA'.
        command (str, optional): [description]. Defaults to 'command'.
        time (int, optional): [description]. Defaults to 4.
        mem (int, optional): [description]. Defaults to 10.
        cpu (int, optional): [description]. Defaults to 8.
        tasks (int, optional): [description]. Defaults to 1.
        dep (str, optional): [description]. Defaults to ''.

    Returns:
        [int]: [Slurm sbatch ID]
    """
    try:
        if dep != '':
            
            dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

        sbatch_command = "sbatch -J {} -o {}.out -e {}.err -t {}:00:00  --mem={}000 --cpus-per-task={} --ntasks={} --wrap='{}' {}".format(
            job_name, sout, serror, time, mem, cpu, tasks,  command, dep)
        
        sbatch_response = subprocess.getoutput(sbatch_command)

        job_id = sbatch_response.split(' ')[-1].strip()

    except Exception:

        log.error("Job submission failed")

    return job_id

def check_status(job_id):
    """
    This function is check status of slurm job

    Args:
        job ([int]): [slurm job id]

    Returns:
        [bool]: [If job completed return true]. Default False.
    """
    d = subprocess.check_output('squeue -j '+str(job_id), shell=True, universal_newlines=True)

    data = list(re.split("\s+ ",d))

    if len(data)==6:

        return True

    return False

def get_cpu():

    """
    This function get actual CPU count of the system 
    Returns:
        rtype: int with 80 % of CPU count
    """

    return math.floor(psutil.cpu_count()*0.8)

def replace_cpu(args, args2):
    '''
    This function replace the actual CPU
    '''

    mat = [i for i in args if any(j in i for j in args2)]
 

    opt , num = mat[0].split(" ")
    count = get_cpu()

    if int(num) > count:
        num = count
        log.warning("number of threads changed to available %s",count)
        
    mat2 = ' '.join([opt,str(num)])
    data = []
    for i in args:
        for j in args2:
            if j in i:
                i= mat2
        data.append(i)

    return data

def change_attribute(args):

    item = ['-g']
    mat = [i for i in item if any(j in i for j in args)]
 

    opt , attr = mat[0].split(" ")
    

    if attr =='ID' :
        attr = 'gene_id'
        log.warning("GTF file provide changing attribute to gene_id")
        
    mat2 = ' '.join([opt,attr])
    data = []
    for i in args:
        for j in args:
            if j in i:
                i= mat2
        data.append(i)

    return data

def get_basename(filePATH):
    """
    This function get the base name of the file from full path 

    Args:
        file ([type]): [path to file]
    """

    return os.path.basename(filePATH)


def get_directory(filePATH):
    """
    This function retrun directory of a file 

    Args:
        file ([type]): [Path to file]

    """

    return os.path.dirname(filePATH)


def get_parent(filePATH):
    """
    This function return the file name without extension

    Args:
        filePATH ([type]): [description]
    """

    return os.path.splitext(filePATH)[0]


def get_file_extension(filePATH):
    """
    This function return the extension of file 

    Args:
        filePATH ([type]): [description]
    """

    return os.path.splitext(filePATH)[1]


def make_directory(dir):
    """
    This function create a directory 

    Args:
        dirName ([type]): [description]
    """
   
    outputdir = os.path.abspath(dir) 

    if os.path.exists(dir):
        parent, base = os.path.split(outputdir)
        counter = 0
        for sibdir in os.listdir(parent):
            if sibdir.startswith(base +'.'):
                ext = sibdir[len(base)+1:]
                if ext.isdecimal():
                    counter = max(counter, int(ext))
        outdir = os.path.join(parent, base+'.'+str(counter+1))

        os.mkdir(outdir)

        log.info(f"Succesfully created directory {outdir}")
        
    else:
        outdir = outputdir
        os.mkdir(outdir)
        log.info(f"Succesfully created directory {outdir}")
    

    return outdir

def check_files(*args):
    """
    This function check if files exist

    Args:
        a list of files to check
    return:
        retrun true only if all files in list exists
    rtype:
        boolean
    """
    for filepath in args:

        if not os.path.isfile(filepath):

            return False

        elif not filepath:

            return False

    return True


def check_path(*args):

    """
    This function check if directory exist

    Args:
        a list of PATH to check
    return:
        retrun true only if all path in list exists
    rtype:
        boolean
    """
    flag = False
    
    for dirPath in args:

        if not (os.path.exists(dirPath) and os.path.isdir(dirPath)):

            flag = True

    if flag:

            return False

    return True

def getFiles(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def findFiles(searchPATH = None, searchPattern = None, recursive = False, verbose = False ):
    """[summary]

    Args:
        searchPATH ([type], optional): [description]. Defaults to None.
        pattern ([type], optional): [description]. Defaults to None.
        recursive (bool, optional): [description]. Defaults to False.
        verbose (bool, optional): [description]. Defaults to False.
    """

    searchResult = []

    if searchPATH is None:

        searchPATH = "./"

    if not check_path(searchPATH):

        return searchResult

    pattern = re.compile(searchPattern)

    if recursive:

        for rootdir, directory, files in os.walk(searchPATH):

            for file in files:

                if bool(pattern.search(file)):

                    searchResult.append(os.path.join(rootdir,file))

        return searchResult
    
    for file in os.listdir(searchPATH):

        if bool(pattern.search(file)):

            searchResult.append(os.path.join(searchPATH,file))
    
    return searchResult

def getGenes(file, combinations, multisheet=True, outDir='pySeqRNA_results'):

    if os.path.exists(outDir):

        make_directory(os.path.join(outDir,"diff_genes"))

    else:

        make_directory("diff_genes")
        outDir = "."

    if multisheet:
        
        for c in combinations:
            df = pd.read_excel(file, sheet_name=c)

            gene = df['Gene'].copy()
            gene = gene.str.replace('gene:','')
            gene.to_csv(os.path.join(outDir,"diff_genes",f"{c}.txt"), sep="\t", index = False)
   
    else:

        for c in combinations:

            df = pd.read(file)

            gene = df['Gene'].copy()

            gene.to_csv(os.path.join(outDir,"diff_genes",f"{c}.txt"), sep="\t", index = False)

    return