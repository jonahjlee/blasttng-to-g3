# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

from data_loader import config, dlib, RoachPass, RoachID, ScanPass

from spt3g import core
import os


if __name__ == '__main__':

    class ScanFrameGenerator(core.G3Module):
        def __init__(self, roach_ids: tuple[int]=None):
            self.roach_ids: tuple[int] = roach_ids if self.roach_ids is not None else (1, 2, 3, 4, 5)
            self.roaches: dict = self._load_roaches()

        def _load_roaches(self):
            """Loads a RoachPass objet for each Roach ID in self.roaches

            Loads all 3 passes by default
            """
            roaches = {roach_id : RoachPass(RoachID(roach_id), ScanPass.ALL, use_rejects_file=False)
                       for roach_id in self.roach_ids}
            return roaches

        def Process(self, frame):
            ...

    out_dir = os.path.join(config.g3_dir, config.version_dir)
    os.makedirs(out_dir, exist_ok=True)

    pipe = core.G3Pipeline()
    pipe.Add(ScanFrameGenerator)
    pipe.Add(core.G3Writer, filename=os.path.join(out_dir, 'testing.g3'))
    pipe.Run()
