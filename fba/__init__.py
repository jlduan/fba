# __init__.py

import subprocess
from pathlib import Path

__version__ = '0.0.12'

if (Path(__file__).resolve().parent.parent / '.git').exists():

    cmd = 'git describe --tags --always --dirty'.split()
    try:
        __version__ = subprocess.check_output(
            cmd, stderr=subprocess.DEVNULL).decode('utf-8').strip()
    except subprocess.CalledProcessError:
        print('Unable to get version number from git tags')
        exit(1)
