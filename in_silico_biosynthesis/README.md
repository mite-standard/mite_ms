In silico biosynthesis example
=========

Organizes code for generating the in silico biosynthesis example display item.

This script demonstrates that the bottromycin A2 precursor peptide can be transformed into the mature product by applying a sequence of reaction rules. 

Reaction rules are retrieved from MITE using the MITE API.

The script produces SMILES of intermediate products, which were used to draw structures using Ketcher, and assembled into a figure.

For more information on MITE, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

## Run the script

*Nota bene*: This installation is only tested on (Ubuntu) Linux.

- Install `python 3.12.x`
- Install hatch (e.g. with `pipx install hatch`)
- Run hatch `hatch env create` to install the dependencies
- Run the script using `hatch run main`

