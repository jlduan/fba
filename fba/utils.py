# utils.py

import logging
import gzip
import bz2
import subprocess
from pathlib import Path


def open_by_suffix(file_name, mode='r'):
    """Opens file based on suffix."""

    # noqa modified from https://stackoverflow.com/questions/18367511/how-do-i-automatically-handle-decompression-when-reading-a-file-in-python
    file_name = str(file_name)

    if file_name.endswith('gz'):
        handle = gzip.open(filename=file_name, mode=mode + 't')
    elif file_name.endswith('bz2'):
        handle = bz2.open(filename=file_name, mode=mode + 't')
    else:
        handle = open(file=file_name, mode=mode)

    return handle


def open_by_magic(file_name):
    """Opens file based on magic."""

    magic_dict = {'\x1f\x8b\x08': (gzip.open, 'rb'),
                  '\x42\x5a\x68': (bz2.BZ2File, 'r')}

    max_len = max(len(x) for x in magic_dict)

    with open(file=file_name) as f:
        file_start = f.read(max_len)

    for magic, (fn, flag) in magic_dict.items():
        if file_start.startswith(magic):
            return fn(file_name, flag)

    return open(file_name, mode='r')


def get_binary_path(binary_name):
    """Gets executable path.

    Parameters
    ----------
    binary_name : str
        The name of the executable.

    Returns
    -------
    str
        The path and name of the executable.

    Raises
    ------
    FileNotFoundError
        If the executable is not found in the PATH.
    """

    binary_path = subprocess.run(
        ['which', binary_name], stdout=subprocess.PIPE, universal_newlines=True
    ).stdout.rstrip()

    if Path(binary_path).is_file():
        return binary_path

    else:
        raise FileNotFoundError(binary_name, 'not found in PATH\n')


def run_executable(cmd_line):
    """Runs executable."""

    proc = subprocess.Popen(
        cmd_line,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    try:
        outs, errs = proc.communicate(timeout=15)
    except subprocess.TimeoutExpired:
        # logger.critical(e, exc_info=True)
        proc.kill()
        outs, errs = proc.communicate()

    return outs, errs


def get_logger(logger_name, log_file=False):
    """Creates a custom logger."""

    FORMATTER = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    logger = logging.getLogger(logger_name)
    logger.setLevel(level=logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(FORMATTER))
    logger.addHandler(console_handler)

    if log_file:
        file_handler = logging.FileHandler(filename='fba.log')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(FORMATTER))
        logger.addHandler(file_handler)

    return logger


def parse_bowtie2_version():
    """Parses bowtie2 version."""

    cmd = [get_binary_path(binary_name='bowtie2'), '--version']
    outs, _ = run_executable(cmd_line=cmd)

    return outs.split(' version ')[1].split()[0]


def parse_samtools_version():
    """Parses samtools version."""

    cmd = [get_binary_path(binary_name='samtools'), '--version']
    outs, _ = run_executable(cmd_line=cmd)

    return outs.split(' ')[1].split()[0]


def parse_kallisto_version():
    """Parses kallisto version."""

    cmd = [get_binary_path(binary_name='kallisto'), 'version']
    outs, _ = run_executable(cmd_line=cmd)

    return outs.rstrip().split(' ')[-1]


def parse_bustools_version():
    """Parses bustools version."""

    cmd = [get_binary_path(binary_name='bustools'), 'version']
    outs, _ = run_executable(cmd_line=cmd)

    return outs.rstrip().split(' ')[-1]
