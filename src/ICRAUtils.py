from functools import wraps
from StringIO import StringIO
from datetime import datetime
from logging import StreamHandler, Filter, WARNING
import logging
from subprocess import check_output as sbcheck_output
import subprocess
import os
from gzip import open as opengz
from os.path import splitext

log_ = logging.getLogger('SGVFinder')
import sys

def timeit(log_, text = None, level = logging.DEBUG):
    def _inner(method):
        '''A decorator which logs the time it took a method to run into debug'''
        @wraps(method)
        def timed(*args, **kw):
            ts = datetime.now()
            result = method(*args, **kw)
            te = datetime.now()
            log_.log(level, ('Method %r took %s' if text is None else text)% (method.__name__, te - ts))
            return result
        return timed
    return _inner

class MaxLevelFilter(Filter):
    '''A filter for the logging module which makes sure only messages below some level pass'''
    def __init__(self, level):
        self.level = level

    def filter(self, record):
        return record.levelno <= self.level
    
def _set_logging_handlers(MIN_LEVEL = logging.INFO, write_to_file=True, warnings_to_stream=True):
    '''Sets handlers printing warning and below to stdout and above warning to stderr'''
    formatter = logging.Formatter('%(asctime)s, %(levelname)s: %(message)s', datefmt='%d/%m %H:%M')
    lower_than_warning = MaxLevelFilter(WARNING)
    rootLogger = logging.getLogger()
    rootLogger.setLevel(MIN_LEVEL)
    
    #Anything that is not an error goes to stdout
    stdout_hdlr = StreamHandler(sys.stdout)
    stdout_hdlr.addFilter(lower_than_warning)
    stdout_hdlr.setFormatter(formatter)
    stdout_hdlr.setLevel(MIN_LEVEL)
    rootLogger.addHandler(stdout_hdlr)
    
    #errors go to stderr
    stderr_hdlr = StreamHandler(sys.stderr)
    stderr_hdlr.setLevel(max(MIN_LEVEL, WARNING) + 1)
    stderr_hdlr.setFormatter(formatter)
    rootLogger.addHandler(stderr_hdlr)
    
    #all logging is also printed to a file
    if write_to_file: 
        log_file_path = ('SGVF_log_%s.txt' % datetime.now()).replace(' ', '_')
        file_hdlr = logging.FileHandler(log_file_path)
        file_hdlr.setFormatter(formatter)
        rootLogger.addHandler(file_hdlr)
    
    #and warnings and above go to a variable
    if warnings_to_stream:
        log_capture_string = StringIO()
        ch = logging.StreamHandler(log_capture_string)
        ch.setLevel(logging.WARNING)
        ch.setFormatter(formatter)
        rootLogger.addHandler(ch)
    
    return log_file_path if write_to_file else None, log_capture_string if warnings_to_stream else None

def _shell_command(com, loglevel = None):
    '''
    Executes a command on the system.
    Will log both the command and its output to the loglevel produced (if not None).
    Differs from os.system and subprocess as both stdout and stderr are returned and are not erroneously printed
    '''
    
    if loglevel is not None:
        log_.log(loglevel, com)
    errnout = ''
    try:
        errnout = sbcheck_output(com, shell = True, stderr = subprocess.STDOUT)
    finally:
        if errnout != '' and loglevel is not None:
            log_.log(loglevel, errnout)
    return errnout

def _tryrm(f):
    try:
        os.remove(f)
    except:pass
    
def _open_gz_indif(filepath):
    return opengz(filepath) if splitext(filepath)[1] in {'.gz', '.gzip'} else open(filepath)

def _none_iter():
    while True:
        yield None