import numpy as np
import g3_utils as ut
import matplotlib.pyplot as plt
from spt3g import core
import pathlib
import time
from astropy.coordinates import SkyCoord

# ============================================================================ #
# CONFIG
# ============================================================================ #

control_computer_g3_dir = pathlib.Path("/media/player1/blast2020fc1/blasttng_g3")
roach1_pass3_file = control_computer_g3_dir / "testing/roach1_pass3.g3"

num_iters = 0

# lambda to get filenames under which to store pipeline outputs at each iteration
iter_file = lambda i: str(control_computer_g3_dir / f"mapmaking/pipeline/stage{i}.g3")

# set exclude_kids to a sequence of KID strings (e.g. ["roach1_0000", ...]) to exclude all kids in the sequence
# set exclude_kids to None or an empty sequence to use all detectors (minus exclude_kids)
kid_rejects = [11, 21, 57, 81, 111, 127, 147, 148, 149, 150, 151, 152, 154, 158, 170, 181, 182, 192, 195,
               211, 223, 225, 245, 263, 282, 286, 293, 297, 319, 327, 331, 333, 334, 336, 337, 340, 341, 349, 352]
exclude_kids = [ut.kid_string(reject_id, roach_id=1) for reject_id in kid_rejects]

# set select_kids to a sequence of KID strings (e.g. ["roach1_0000", ...]) to exclude all kids not in the sequence
# if select_kids is None, no additional filtering is applied
select_kids = None

source = SkyCoord.from_name('RCW92')
map_args = {
    'ra0': source.icrs.ra.deg * core.G3Units.deg,
    'dec0': source.icrs.dec.deg * core.G3Units.deg,
    'xlen': 2.5 * core.G3Units.deg,
    'ylen': 1.3 * core.G3Units.deg,
    'res': 1 * core.G3Units.arcmin,
}

if __name__ == "__main__":

    # ============================================================================ #
    # PREPROCESSING STAGE
    # ============================================================================ #

    print("\nBeginning Pre-Processing Stage...\n")
    start_time = time.perf_counter()

    stats = ut.DetectorStats(data_key="df")

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=str(roach1_pass3_file))
    # create non-normalized df & save
    pipe.Add(ut.add_cal_lamp_df, iq_key="cal_lamp_data", exclude_kids=exclude_kids, select_kids=select_kids)
    pipe.Add(ut.AddScanDF, exclude_kids=exclude_kids, select_kids=select_kids)
    pipe.Add(core.G3Writer, filename=iter_file(0))  # note: iteration 0 file does not have ra/dec or normalized df
    # compute medians & standard deviation for each detector scan
    pipe.Add(stats)
    pipe.Run(profile=True)

    detector_medians = np.median(np.array(stats.medians), axis=0)
    detector_stds = np.median(np.array(stats.stds), axis=0)

    # ============================================================================ #
    # INITIAL MAP
    # ============================================================================ #

    print("\nBeginning Map-Making Stage (Naive Common-Mode Removal)...\n")

    binner0 = ut.MapBinner(timestreams="it1", stds=detector_stds, **map_args)

    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=iter_file(0))
    pipe.Add(ut.add_radec_so3g)
    pipe.Add(ut.NormalizeDF, detector_medians=detector_medians)
    pipe.Add(ut.remove_common_mode, in_key="norm_df", out_key="it1")
    pipe.Add(core.G3Writer, filename=iter_file(1))
    pipe.Add(binner0)
    pipe.Run(profile=True)

    # ============================================================================ #
    # COMMON-MODE ITERATION
    # ============================================================================ #

    prev_binner = binner0
    prev_iter = 1

    for _ in range(num_iters):
        print(f"\nIterating Common-Mode Removal...\n")

        with np.errstate(invalid='ignore'):
            prev_map = prev_binner.data / prev_binner.hits

        iter_binner = ut.MapBinner(timestreams=f"it{prev_iter}", stds=detector_stds, **map_args)

        pipe = core.G3Pipeline()
        pipe.Add(core.G3Reader, filename=iter_file(prev_iter))
        pipe.Add(ut.common_mode_iter,
                 kid_shifts=prev_binner.kid_shifts,
                 in_key=f"it{prev_iter}", out_key=f"it{prev_iter + 1}",
                 prev_map=prev_map, **map_args)
        pipe.Add(iter_binner)
        pipe.Add(core.G3Writer, filename=iter_file(prev_iter + 1))
        pipe.Run(profile=True)

        prev_binner = iter_binner
        prev_iter += 1

    end_time = time.perf_counter()

    # ============================================================================ #
    # OUTPUT
    # ============================================================================ #

    # determine the duration of the scan data
    ffg = ut.FirstFrameGrabber(core.G3FrameType.Scan)
    lfg = ut.LastFrameGrabber(core.G3FrameType.Scan)
    pipe = core.G3Pipeline()
    pipe.Add(core.G3Reader, filename=str(roach1_pass3_file))
    pipe.Add(ffg)
    pipe.Add(lfg)
    pipe.Run()
    scan_duration = (lfg.last_frame["az"].stop.time - ffg.first_frame["az"].start.time) / core.G3Units.s

    # print out some stats for the map making run
    print("\nCreated a map with:")
    print(f"    {len(prev_binner.kids)} detectors")
    print(f"    {scan_duration}s of data.")
    print(f"    {num_iters} commmon-mode iterations (+1 naive iteration at the start)")
    print(f"Total Time: {end_time - start_time}s")

    # save the map (comment this out if this is unwanted)
    prev_binner.plot(show=False)
    plt.savefig(pathlib.Path(__file__).parent / 'map.png')
