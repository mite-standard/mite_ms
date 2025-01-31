MITE sequence similarity network
=========

Organizes code for generating the files for a sequence similarity network of MITE entries.

The script prepares files for a specific record of mite_data (a version), given as command line argument.

The generated fasta file still needs to be run with EFI-EST and visualized in Cytoscape. The generated metadata.csv file can be used to annotate this file.

For more information on MITE, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

## Run the script

*Nota bene*: This installation is only tested on (Ubuntu) Linux.

- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Run hatch `hatch env create` to install the dependencies
- Run the script using `hatch run main` to see the command line options

