MITE tmap reaction SMILES
=========

Organizes code for generating a tmap plot of reaction diversity of MITE

The script prepares files for a specific record of mite_data (a version), given as command line argument.

The script converts MITE example reaction substrate->products into reaction SMILES, calculates drfp fingerprints, and visualizes them with a tmap.

For more information on MITE, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

## Run the script

*Nota bene*: This installation is only tested on (Ubuntu) Linux.

- Install `python 3.12.x`
- Install (micro)mamba
- Run the following commands to install the dependencies
```commandline
micromamba create -n mite_tmap python=3.9
micromamba activate mite_tmap
micromamba install -c tmap tmap
micromamba install conda-forge::python-annoy
pip install -r requirements.txt
```
- Run the script using `python3 module/main.py` to see the command line options
- Remove the environment:
```commandline
micromamba deactivate
micromamba env remove -n mite_tmap
```
