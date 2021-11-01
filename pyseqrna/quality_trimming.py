#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Title: This modules contains read quality check functions for pySeqRNA
Created : 
@author : Naveen Duhan
'''

import os
import shutil
import subprocess
import pkg_resources
from pyseqrna import PyseqrnaLogger
from pyseqrna import pyseqrna_utils as pu

# Intialize the logger

pl = PyseqrnaLogger(mode='a')