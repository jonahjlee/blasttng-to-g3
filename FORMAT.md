# Data Format

This document describes the format of data in the produced BLAST-TNG ``.g3`` files.
In particular, descriptions are provided for:
- The folder structure of outputted ``.g3`` files
- Which G3Frames are in present and how often
- Which keys are present in each frame, and details about the contents.

Note: The data format is a work in progress and is subject to change

## Folder Structure
TBD

## Frames

### 1. Scan Frames 

The [CCAT rfsoc-streamer README](https://github.com/ccatobs/rfsoc-streamer/blob/main/README.md)
lists the keys present in CCAT scan frames. The first ``.g3`` files I make will contain a subset of these keys, 
sufficient for our mapmaking purposes.

### 2. Observation Frames

TBD

### 3. Calibration Frames

TBD

### 4. Misc
- [PipelineInfo](https://cmb-s4.github.io/spt3g_software/frames.html#pipelineinfo) Frames
  - PipelineInfo frames are added automatically by G3Pipeline.

- [EndProcessing](https://cmb-s4.github.io/spt3g_software/frames.html#endprocessing) Frames
  - EndProcessing frames are added automatically by G3Pipeline.

## Keys / Frame Objects
TBD