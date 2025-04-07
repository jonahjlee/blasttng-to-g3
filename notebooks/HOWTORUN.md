# Remote Jupyter Notebooks

## Purpose

Since G3 files will mostly be stored on the CCAT Control Computer (in ``/media/player1/blast2020fc1/blasttng_g3``),
it is useful to run notebooks in a Jupyter server on the control computer, then to use SSH tunneling to interact with
the Jupyter notebook through the client's browser. This minimizes the tedious process of saving and transferring data that is needed
for anything other than text output (e.g. Matplotlib plots).

## How To Run

Run these commands on the desired client machine.
For this walkthrough, it assumed that the client is already on the CCAT [tailnet](https://tailscale.com/).

1. SSH into the control computer. Follow instructions to authenticate if needed.
    ```
    ssh player1@cube
    ```
2. Enter the project repository on the remote machine
    ```
    cd ~/map_making_jonah/blasttng-to-g3
    ```
3. Set up ``PYTHONPATH`` to fix possible import issues (optional)
    ```
    export PYTHONPATH=$PYTHONPATH:$PWD
    ```
4. Activate the virtual environment
    ```
    source .venv/bin/activate
    ```
5. Start the Jupyter server as a module. Running it this way will ensure the virtual environment is used rather than the global python installation.
    ```
    python3 -m jupyter notebook --no-browser --port 8888
    ```
6. Next, copy the URL that appears in the output (e.g. http://localhost:8888/tree?token=placeholder123456). Return to the terminal on your client machine for the following commands. It may be useful to do this in a new terminal. This will forward the port 8888 on the client to the localhost of the __remote__ control computer at port 8888.
    ```
    # on local machine (not control computer)
    ssh player1@cube -NL 8888:localhost:8888
    ```
7. Finally paste the link copied earlier into your browser and you should be good to go!