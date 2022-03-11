# __init__.py

import subprocess
from pathlib import Path

__version__ = '0.0.12'

package_root = Path(__file__).resolve().parent.parent
if (package_root / '.git').exists():

    cmd = 'git describe --tags --always --dirty'.split()
    try:
        __version__ = subprocess.check_output(
            cmd, cwd=package_root).decode().strip()

    except subprocess.CalledProcessError:
        print('Unable to get version number from git tags')
        exit(1)
