# BLAST-TNG to G3
Tools for re-packaging BLAST-TNG data into g3 files.

## Context

A Quick-Look Mapmaker is in development for CCAT, with the purpose of quickly producing maps and enabling time-domain science. Proof-of-concept code can be found at:
- The main [map_making_blasttng](https://github.com/ccat850/map_making_blasttng) repository, which contains the ``map_maker_iterative`` project
- Jonah's [map_making_blasttng_jl](https://github.com/jonahjlee/map_making_blasttng_jl) fork, which extends ``map_maker_iterative`` by allowing for selection of multiple roaches and specific passes of RCW-92.

In its current state, ``map_maker_iterative`` accesses the (large volume of) BLAST-TNG RCW-92 observation data by loading ``.npy`` files from the CCAT Control Computer's local file system. These numpy files themselves are the result of processing of dirfiles using the code in https://github.com/ccat850/map_making_blasttng/blob/main/dirfile_to_npy.py. _However_, the actual detector data for CCAT will be in the form of .g3 files as outputted by https://github.com/ccatobs/rfsoc-streamer. This means that as the Quick-Look Mapmaker is re-written into its production form, it must be built to work with .g3 files -- and since we do not currently have complete sample data in the form we expect for CCAT, sample .g3 files from BLAST-TNG data will be useful to enable development and testing of the new mapmaker. Below, the motivation for this is outlined in more detail.

## Motivation

- Detector data for CCAT will be available to the Quick-Look Mapmaker in the form of .g3 files
  - Packaging BLAST-TNG data into .g3 files provides ideal sample data:
  - We already know how to make maps with BLAST-TNG data ([main repo](https://github.com/ccat850/map_making_blasttng) / [Jonah’s fork](https://github.com/jonahjlee/map_making_blasttng_jl))
  - The data shares some similarity with CCAT:
    - MKID data on several RF networks
    - Relatively large number of detectors (though still much lower than CCAT!)
  - BLAST-TNG data is already available & ready to go on the CCAT Control Computer
- Once packaged, the BLAST-TNG .g3 can be used for development and testing of a new quick-look mapmaker
  - We can begin porting our mapmaking code to the .g3 format
  - We can experiment with pipelines and various ways to improve performance
- We can also create synthetic .g3 files with more data to determine performance characteristics of the mapmaker (i.e., how does time/space scale with number of detectors or downsampling rate?)

## Goals

- We want this sample .g3 data to be as close as possible to the format we will get for CCAT
  - The closer it is, the lower the cost will be once we transition from BLAST-TNG .g3 data to real CCAT data.
  - Communicating with DAQ team will help establish parity between the formats
- It would be helpful to have several different forms/versions of .g3 sample data
  - Various data rates/file sizes due to different number of detectors
  - Could duplicate detectors to create more data, or just add garbage data for performance testing purposes
    - duplicating data may lead to inaccurate performance results if some optimizations cause unnatural speedups
- Provide adequate documentation (slides, READMEs, code comments, docstrings…) such that somebody in the future can:
  - Create new .g3 files with different data
  - Understand the format of the existing BLAST-TNG .g3 files and work with them

## Project Structure

The ``g3_packager`` module is responsible for generating the G3 files, and is intended to be run on the CCAT Control computer.
Running ``main.py`` will read the data required using the ``data_loader`` package and create a new G3 file for the data.

The ``data_loader`` package consists of parts of code adapted from [map_making_blasttng_jl](https://github.com/jonahjlee/map_making_blasttng_jl) that are required
to load (and minimally process) BLAST-TNG data.

This project also contains a map-maker that effectively ports the functionality from [map_making_blasttng_jl](https://github.com/jonahjlee/map_making_blasttng_jl) to
make maps using G3 pipelines. The main script is in ``profiling/mapmaker.py``. The output from ``mapmaker.py`` is given
by the quotient ``binner.data / binner.hits`` for the binner which took in the  most recently processed data.
