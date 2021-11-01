#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: naveen duhan

Read pyrpipe configuration
"""
import os 
import logging

def PyseqrnaLogger(mode):

    """[intialize logger in the pySeqRNA modules]

    Args:
        logger (str): [logger name for the module]
        logFile (str): [file name for logging]
    """
    logger = logging.getLogger(__name__)
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
