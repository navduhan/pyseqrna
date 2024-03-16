#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: This modules contains utility functions for pySeqRNA

:Created : July 11, 2021

:Author : Naveen Duhan
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
from argparse import Namespace
from openpyxl import load_workbook

def PyseqrnaLogger(mode, log):

    """ This function intialize logger in the pySeqRNA modules

    :param  mode: Logger name for the module

    :param  log: File name for logging
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

def read_runconfig(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    # Create Namespace object to hold the values
    options = Namespace()

    # General section
    general = config['General']
    options.input_file = general.get('input_file', '')
    options.samples_path = general.get('samples_path', '')
    options.reference_genome = general.get('reference_genome', '')
    options.feature_file = general.get('feature_file', '')
    options.version = general.get('version', '')

    # Internal section
    internal = config['Internal']
    options.source = internal.get('source', 'ENSEMBL')
    options.species = internal.get('species', '')
    options.speciestype = internal.get('organismType', 'plants')
    options.outdir = internal.get('outdir', 'pySeqRNA_results')
    options.paired = internal.getboolean('paired', fallback=False)
    options.fastqc = internal.getboolean('fastqc', fallback=False)
    options.fastqc2 = internal.getboolean('fastqc2', fallback=False)
    options.ribosomal = internal.getboolean('ribosomal', fallback=False)
    options.rnadb = internal.getboolean('rnadb', fallback=False)
    options.mmgg = internal.getboolean('multimappedGroups', fallback=False)
    options.minreadcounts = int(internal.get('minReadCounts', 100))
    options.percentsample = float(internal.get('percentSample', 0.5))
    combination_value = internal.get('combination', 'all')
    options.combination = combination_value.split() if combination_value != 'all' else combination_value
    options.fdr = float(internal.get('fdr', 0.05))
    options.fold = int(internal.get('fold', 2))
    options.noreplicate = internal.getboolean('noreplicate', fallback=False)
    options.normalizecount = internal.get('normalizeCount', 'RPKM')
    options.coexpression = internal.get('coexpression', fallback=False)
    options.heatmap = internal.getboolean('heatmap', fallback=False)
    options.heatmaptype = internal.get('heatmapType', 'counts')
    options.pcaplot = internal.getboolean('pcaplot', fallback=False)
    options.pcaplottype = internal.get('pcaplotType', 'RPKM')
    options.maplot = internal.getboolean('maPlot', fallback=False)
    options.volcanoplot = internal.getboolean('volcanoPlot', fallback=False)
    options.vennplot = internal.getboolean('vennPlot', fallback=False)
    venncombination_value = internal.get('vennCombinations', 'random')
    options.venncombination = venncombination_value.split() if venncombination_value != 'random' else venncombination_value
    options.cluster = internal.getboolean('cluster', fallback=True)

    # FunctionalAnnotation section
    annotation = config['FunctionalAnnotation']
    options.geneontology = annotation.getboolean('geneOntology', fallback=False)
    options.keggpathway = annotation.getboolean('keggPathway', fallback=False)

    # ExternalTools section
    tools = config['ExternalTools']
    options.trimming = tools.get('trimming', 'trim_galore')
    options.aligner = tools.get('aligner', 'STAR')
    options.quantification = tools.get('quantTool', 'featureCounts')
    options.detool = tools.get('deTool', 'DESeq2')

    # Computational section
    computational = config['Computational']
    options.slurm = computational.getboolean('slurm', fallback=False)
    options.slurm_partition = computational.get('slurm_partition', 'compute')
    options.threads = computational.get('threads', get_cpu())
    options.memory = int(computational.get('memory', 16))
    options.param = computational.get('param', fallback=None)
    options.resume = computational.get('resume', 'all')

    return options

def read_input_file(infile, inpath, paired = False):

    """This function reads input sample file and convert into a dictionary. It also make all possible combination for DEG analysis. Target dataframe for differential analysis.

    :param inputFile: input sample file containing the infromation about project

    :param inputPath: Path for input fastq files

    :param pairedEND: Check if reads are paired end]. Defaults to False

    :return:   samples, combinations and targets for differential expression
    :rtype: A dictionary
    """
    samples = {} # dictionary to collect sample information from input sample file
    factors = [] # list for collecting the Identifier infromation from input sample file
    combinations = [] # list for genrating combinations for differential expression analysis

    try:
      
        log.info(f"Reading input file {infile}")

        if infile.endswith('.xlsx'):
            df = pd.read_excel(infile, comment="#")
        elif infile.endswith('.csv'):
            df = pd.read_csv(infile,comment="#")
        elif infile.endswith('.txt'):
            df = pd.read_csv(infile, sep="\s+" , comment="#")
        else:
            raise TypeError("Unsupported input format")

        rows = df.values.tolist()

        for lines in rows:
            

            if paired:
                samples[lines[1]] = [lines[1],lines[2],os.path.join(inpath,lines[3]), os.path.join(inpath,lines[4])]
            else:
                samples[lines[1]] = [lines[1], lines[2], os.path.join(inpath,lines[3])]

            if lines[2] not in factors:
                factors.append(lines[2])


        log.info("Input file %s read succesfully", infile)

    except Exception:

        log.error("Please provide a valid input file")

        sys.exit()

    try:

        # create combinations from factors 
        for i in factors:
            for j in factors:

                if i != j:

                    if j+"-"+i not in combinations:

                        combinations.append(i+"-"+j)
        log.info("Combination created succesfully from %s", infile)
    except Exception:

        log.error("Please provide a valid input file")
    
    samplename = []
    sample = []

    try:
        for k, s in samples.items():

            samplename.append(k)

            sample.append(s[1])

        targets = pd.DataFrame(samplename,index=[i for i in samplename])

        targets = targets.assign(sample=sample)

        log.info("targets dataframe for differenatial created succesfully from %s", infile)

    except Exception:

        log.error("Please provide a valid input file")

    return {'samples': samples, "combinations": combinations, "targets": targets}



def parse_config_file(infile):
    """
    This function parse the config file for all the programs used in pySeqRNA

    :param: configFile: <program>.ini config file containing arguments. 

    :retrun: Program specific arguments 

    :rtype: a dictionary
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

        log.info("Config generated succesfully from %s", infile)

    except Exception:

        log.error("Please provide a valid config file")

    return sections_dict

def clusterRun(job_name='pyseqRNA',sout=" pyseqrna", serror="pyseqrna", partition='compute', command='command', time=4, mem=10, cpu=8, tasks=1, dep=''):
    """
    This function is for submitting job on cluster with SLURM job Scheduling

    :param job_name: Slurm job name on HPC. Defaults to 'pyseqRNA'.

    :param command:  Command to excute on HPC.

    :param time: Slurm Job time allotment. 

    :param mem: Memory to use in GB.

    :param cpu: Number of CPUs to use for the job.

    :param tasks: Number of tasks to execute.

    :param dep: Slurm Job dependency. Defaults to ''.

    :returns:

        :rtype: Slurm sbatch ID
    """
    try:
        if dep != '':
            
            dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

        sbatch_command = "sbatch -J {} -o {}.out -e {}.err -t {}:00:00 -p {} --mem={}000 --cpus-per-task={} --ntasks={} --wrap='{}' {}".format(
            job_name, sout, serror, time, partition, mem, cpu, tasks,  command, dep)
        
        sbatch_response = subprocess.getoutput(sbatch_command)

        job_id = sbatch_response.split(' ')[-1].strip()

    except Exception:

        log.error("Job submission failed")

    return job_id

def check_status(job_id):
    """
    This function is check status of slurm job

    :param job_id: slurm job id

    :returns: True/False

    :rtype: If job completed return True. Default False.
    """
    d = subprocess.check_output('squeue -j '+str(job_id), shell=True, universal_newlines=True)

    data = list(re.split("\s+ ",d))

    if len(data)==6:

        return True

    return False

def get_cpu():

    """
    This function get actual CPU count of the system 

    :returns: Integer 
    
    :rtype: int with 80 % of CPU count
    """

    return math.floor(psutil.cpu_count()*0.8)

def replace_cpu(args, args2):
    '''
    This function replace the actual CPU in config file.

    :returns: Change CPU count to 80% of available CPU 
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

    '''
    This function changes the attribute for feature counts based on GFF or GTF file.
    '''

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

    :param  filePATH: Path to file.
    """

    return os.path.basename(filePATH)


def get_directory(filePATH):
    """
    This function retrun directory of a file 

    :param  filePATH: Path to file.
    """

    return os.path.dirname(filePATH)


def get_parent(filePATH):
    """
    This function return the file name without extension

    :param  filePATH: Path to file.
    """

    return os.path.splitext(filePATH)[0]


def get_file_extension(filePATH):
    """
    This function return the extension of file 

    :param  filePATH: Path to file.
    """

    return os.path.splitext(filePATH)[1]


def make_directory(dir, dryrun=False):
    """
    This function create a directory 

    :param  dir: Directory name.

    :returns: Name of created directory.
    """
    if dryrun:

        return dir
    
    else:
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

    :param args: List of files to check.

    :return: True or False.
 
    :rtype: Retrun true only if all files in list exists in a directory.
    
    """
    for filepath in args:

        if not os.path.isfile(filepath):

            return False

        elif not filepath:

            return False

    return True


def check_path(*args):

    '''
    This function check if directory exist

    :param args: List of directory to check.

    :return: True or False
 
    :rtype: Retrun true only if all files in list exists in a directory.
    '''
    flag = False
    
    for dirPath in args:

        if not (os.path.exists(dirPath) and os.path.isdir(dirPath)):

            flag = True

    if flag:

            return False

    return True

def getFiles(pattern, path):
    """
    This function searches all files containing patterns in a directory.

    :param pattern: Patteren to search.

    :param path: A direcory path.

    :returns: All files containing a pattern.
    """
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def findFiles(searchPATH = None, searchPattern = None, recursive = False, verbose = False ):
    """
    This function find searches files.

    :param searchPATH: Search directory.

    :param pattern: Pattern to search.

    :param recursive: True if want to search recursively. Defaults to False.

    :param verbose: True if want to print output. Defaults to False.
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

def getGenes(file=None, combinations=None, multisheet=True,  geneType='all',  outDir='.', mmg=False):

    """
    This function extract genes from filtered differentiall expressed genes.

    :param file: Filtered DEGs file.

    :param combinations: Comparison list contaning samples.

    :param multisheet: True if file is multisheet.

    :geneType: Genes to extract all, up , down. Default is all.

    :param outDir: Output directory. Default is current working directory.
    """

    if os.path.exists(outDir):

        if os.path.exists(os.path.join(outDir,"diff_genes")):
            out = os.path.join(outDir,"diff_genes")
            
        else:
            out = make_directory(os.path.join(outDir,"diff_genes"))

    else:

        if os.path.exists(os.path.join(outDir,"diff_genes")):
            out = os.path.join(outDir,"diff_genes")
        else:
            out = make_directory(os.path.join(outDir,"diff_genes"))


    if multisheet:
        
        for c in combinations:


            df = pd.read_excel(file, sheet_name=c)
            if mmg:
                mmg_genes = df['Gene'].values.tolist()

                mgenes = []
                for m in mmg_genes:
                    mp = m.split("-")
                    for p in mp:
                        mgenes.append(p)
                mgenesDF = pd.DataFrame(mgenes, columns=['Gene'])
                if geneType == 'all':
                    mgenesDF.to_csv(os.path.join(out,f"{c}_mmg.txt"), sep="\t", index = False)
                if geneType == 'up':
                    mgenesDF.to_csv(os.path.join(out,f"{c}_mmg_up.txt"), sep="\t", index = False)
                if geneType == 'down':
                    mgenesDF.to_csv(os.path.join(out,f"{c}_mmg_down.txt"), sep="\t", index = False)
            else:
                gene = df['Gene'].copy()
                gene = gene.str.replace('gene:','').str.upper()
                if geneType == 'all':
                    gene.to_csv(os.path.join(out,f"{c}.txt"), sep="\t", index = False)
                if geneType == 'up':
                    gene.to_csv(os.path.join(out,f"{c}_up.txt"), sep="\t", index = False)
                if geneType == 'down':
                    gene.to_csv(os.path.join(out,f"{c}_down.txt"), sep="\t", index = False)
   
    else:

        for c in combinations:

            df = pd.read(file)

            gene = df['Gene'].copy()

            gene.to_csv(os.path.join(out,f"{c}.txt"), sep="\t", index = False)

    return  out



def parse_gff(file):
    """
    This function parse a gene feature file in a dataframe for gene IDs

    :param file: A gene feature file.
    """

    gtf_file = pd.read_csv(file ,sep="\t", header=None, comment="#")

    gtf_file.columns = ['seqname', 'source', 'feature', 'start', 'end', 's1','strand', 's2', 'identifier']
    gene = pd.DataFrame(gtf_file[gtf_file['feature'] == 'gene'])

    gene['Gene'] = list(map(lambda x: re.search(r'ID=(.*?);',x,re.MULTILINE).group(1).split("gene-")[1],gene['identifier'].values.tolist()))
    gene['entrez']= list(map(lambda x: re.search(r'GeneID:(.*?);',x,re.MULTILINE).group(1).split(",")[0],gene['identifier'].values.tolist()))
    gene_list = gene[['Gene', 'entrez']]

    # cds = pd.DataFrame(gtf_file[gtf_file['feature'] == 'CDS'])
    # cds['cds'] = list(map(lambda x: re.search(r'ID=(.*?);',x,re.MULTILINE).group(1).split("cds-")[1],cds['identifier'].values.tolist()))
    # cds['entrez']= list(map(lambda x: re.search(r'GeneID:(.*?);',x,re.MULTILINE).group(1).split(",")[0],cds['identifier'].values.tolist()))
    # cds_list = cds[['cds', 'entrez']]

    # final = gene_list.merge(cds_list, on='entrez')

    final = gene_list.drop_duplicates()

    return final

def change_ids(df, file):
    """
    This function changes ids to other ids.

    :param df: A DataFrame containing IDs and synonym IDs.

    :param file: A DataFrame in which IDs needs to be replaced.
    """
    
    pl = df.values.tolist()
    pp = {}
    for p in pl:
        pp[p[1]] = p[0]

    for i, row in file.iterrows():
        Genes = str(file.at[i, 'Genes']).split(",")
        result = []
        for gene in Genes:
            result.append(pp[gene])
        res = ",".join(result)
        file.at[i, 'Genes'] = res
    return file

def add_MMG(degDF=None, anotDF=None, combination=None):

    filtered_Genes = pd.read_excel(degDF, sheet_name=combination)
    Fgenes = filtered_Genes[['MMG', 'Gene']].values.tolist()
    
    Fgenes_List = []

    for k in Fgenes:
        for m in k[1].split("-"):
            Fgenes_List.append([k[0], m.upper()])

    genes = anotDF[['Genes']].values.tolist()
    MMG_genes = []
    for g in genes:
        gg = g[0].split(",")
        p = []
        for g in gg:
            
            for m in Fgenes_List:
                if g == m[1]:
                    p.append(m[0])
        MMG_genes.append(','.join(p))

    
    anotDF['MMG'] = MMG_genes



    return anotDF