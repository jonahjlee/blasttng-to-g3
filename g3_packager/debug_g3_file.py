# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

from data_loader import config
from spt3g import core
from typing import Union

import argparse
import os

from data_loader.config import version_dir

def load(path) -> list[core.G3Frame]:
    return list(core.G3File(path))

if __name__ == '__main__':

    breakpoint()

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path',
                        default=None,
                        description=f"File path to folder containing .g3 file from g3_dir ({config.g3_dir})",)
    args = parser.parse_args()

    if args.filepath is None:
        filepath: os.path = os.path.join(r"C:\Users\Jonah\PycharmProjects\blasttng-to-g3\tmp", version_dir)
    else:
        filepath: os.path = os.path.join(r"C:\Users\Jonah\PycharmProjects\blasttng-to-g3\tmp", args.filepath)

    files = [file for file in os.listdir(filepath) if file.endswith('.g3')]

    print(f"G3 files in {filepath}", files)

    breakpoint()
