# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

import argparse
from fileinput import filename

import so3g
from spt3g import core
import os

class FrameCounter(core.G3Module):
    def __init__(self):
        super(FrameCounter, self).__init__()
        self.previous_type = None
        self.num_repeats = 0
    def Process(self, frame):
        type = frame.type
        if type == self.previous_type:
            self.num_repeats += 1
            print(f"{type} (x{self.num_repeats + 1})", end='\r')
        else:
            print()
            print(f"{type}", end='\r')
            self.num_repeats = 0
        self.previous_type = type

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default=None)
    args = parser.parse_args()

    filename = os.path.normpath(r"/mnt/c/Users/Jonah/PycharmProjects/blasttng-to-g3/tmp/testing/testfile.g3")
    if args.file is not None:
        # override default file
        filename = os.path.normpath(args.file)

    print(f"G3 file: {filename}")

    def load() -> list[core.G3Frame]:
        return list(core.G3File(filename))

    def dump(profile=False) -> None:
        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader, filename=filename)
        pipe.Add(core.Dump)
        pipe.Run(profile=profile)

    breakpoint()
