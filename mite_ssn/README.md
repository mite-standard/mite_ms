MITE sequence similarity network
=========

Organizes code for generating the files for a sequence similarity network of MITE entries.

The script prepares files for a specific record of mite_data (a version), given as command line argument.

The generated fasta file still needs to be run with EFI-EST and visualized in Cytoscape. The generated metadata.csv file can be used to annotate this file.

Further, MITE protein sequences are blasted against NCBI non-redundant protein database and the number of matches per entry determined.
A boxplot is generated to indicate the distributions of BLAST matches against the NCBI NR protein database. 
The statistical summary can be found in an accompanying json file.

For more information on MITE, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

## Run the script

*Nota bene*: This installation is only tested on (Ubuntu) Linux.

- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Run hatch `hatch env create` to install the dependencies
- Run the script using `hatch run main` to see the command line options

