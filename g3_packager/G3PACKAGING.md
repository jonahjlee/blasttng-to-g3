# G3 Packaging

The ``g3_packager`` module is responsible for generating the G3 files, and is intended to be run on the CCAT Control computer.
Running ``main.py`` will read the data required using the ``data_loader`` package and create a new G3 file for the data.

The script reads data from the numpy files stored on the control computer (in ``/media/player1/blast2020fc1/fc1/``),
and also loads detector shifts from ``./detector_layouts``. Currently, ``roach1_shifts_radec.npy`` is the only file
that has detector shifts. To generate more, see ``../notebooks/map_making.ipynb`` to see how source coordinates are found
in SingleMapBinners, processed into angular ``shifts_radec`` and saved.