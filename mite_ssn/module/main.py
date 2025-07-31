"""Prepare files and metdata for MITE sequence similarity network

MIT License

Copyright (c) 2024 to present Mitja M. Zdouc and individual contributors.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import json
import logging
import os
from pathlib import Path
import shutil
import sys
import time

import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from pydantic import BaseModel
import requests
import seaborn as sns

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)


class AbstractManager(BaseModel):
    """Class to contain general attributes

    data: path to mite_data json files
    fasta: path to fasta files
    ncbi_results: path to NCBI NR BLAST results
    output: resulting files for SSN
    """

    data: Path = Path(__file__).parent.joinpath("data/data")
    fasta: Path = Path(__file__).parent.joinpath("data/fasta")
    ncbi_results: Path = Path(__file__).parent.joinpath("ncbi_nr")
    output: Path = Path(__file__).parent.joinpath("output")


class DownloadManager(AbstractManager):
    """Download files from Zenodo record and unpack

    Attributes:
        record: the record to download
        location: the location to download data to
        record: path to the record file
        record_unzip: path to unzipped record file
        version: path to file containing the version of mite_data used
    """

    record: str
    location: Path = Path(__file__).parent.joinpath("data")
    record: Path = Path(__file__).parent.joinpath("data/record.zip")
    record_unzip: Path = Path(__file__).parent.joinpath("data/record")
    version: Path = Path(__file__).parent.joinpath("version.json")

    def run(self) -> None:
        """Call methods for downloading and moving data"""

        if self.location.exists():
            logger.warning(
                "MITE data folder already present - skip download. Remove the folder manually if a new version needs to be downloaded."
            )
            return

        self.location.mkdir(parents=True)
        self.download_data()
        self.organize_data()

    def download_data(self) -> None:
        """Download data from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        response_metadata = requests.get(
            f"https://zenodo.org/api/records/{self.record}"
        )
        if response_metadata.status_code != 200:
            logger.fatal(
                f"Error fetching 'mite_data' record metadata: {response_metadata.status_code}"
            )
            raise RuntimeError

        record_metadata = response_metadata.json()
        version = record_metadata["metadata"]["version"]
        files_url = record_metadata["files"][0]["links"]["self"]

        response_data = requests.get(files_url)
        if response_data.status_code != 200:
            logger.fatal(
                f"Error downloading 'mite_data' record: {response_data.status_code}"
            )
            raise RuntimeError

        with open(self.version, "w") as f:
            f.write(json.dumps({"version_mite_data_used": f"{version}"}))

        with open(self.record, "wb") as f:
            f.write(response_data.content)

    def organize_data(self) -> None:
        """Unpacks data, moves to convenient location, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        shutil.unpack_archive(
            filename=self.record, extract_dir=self.record_unzip, format="zip"
        )
        if not self.record_unzip.exists():
            logger.fatal(f"Could not find the unzipped directory {self.record_unzip}.")
            raise NotADirectoryError

        matching_dirs = list(self.record_unzip.glob("mite-standard-mite_data-*"))
        if not matching_dirs:
            logger.fatal(
                f"Could not determine data storage location in downloaded directory."
            )
            raise RuntimeError

        subdir = matching_dirs[0]

        shutil.move(
            src=self.record_unzip.joinpath(subdir).joinpath("mite_data/data").resolve(),
            dst=self.location.resolve(),
        )

        shutil.move(
            src=self.record_unzip.joinpath(subdir)
            .joinpath("mite_data/fasta")
            .resolve(),
            dst=self.location.resolve(),
        )

        os.remove(self.record)
        shutil.rmtree(self.record_unzip)


class BlastManager(AbstractManager):
    """Class to annotate MITE entries against NCBI-NR using their BLAST API"""

    def run(self):
        self.ncbi_results.mkdir(exist_ok=True)

        for fasta in self.fasta.iterdir():
            if self.ncbi_results.joinpath(f"{fasta.stem}.xml").exists():
                continue

            record = SeqIO.read(fasta, "fasta")
            query = f""">{record.id}
            {record.seq}
            """

            logger.info(f"Submitting BLASTp entry of {record.id}")
            result_handle = NCBIWWW.qblast(
                program="blastp",
                database="nr",
                sequence=query,
                expect=1e-5,
                hitlist_size=5000,  # Maximum results
                format_type="XML",
            )
            raw_blast_output = result_handle.read()
            with open(self.ncbi_results.joinpath(f"{record.id}.xml"), "w") as file:
                file.write(raw_blast_output)

            logger.info("BLAST search completed. Results saved.")

            time.sleep(10)  # limits rate to prevent IP ban by NCBI


class MetadataManager(AbstractManager):
    """Class to generate metadata files for EFI-EST

    metadata_efi_est: a dict containing efi-est formatted metadata
    fasta_efi_est: buffer for concatenating
    """

    metadata_efi_est: dict = {
        "key": [],
        "mite_acc": [],
        "enzyme_name": [],
        "enzyme_description": [],
        "tailoring": [],
        "id_uniprot": [],
        "id_genpept": [],
        "id_mibig": [],
        "ncbi_nr_matches": [],
    }
    fasta_efi_est: list = []
    nr_blast_matches: dict = {
        "mite_acc": [],
        "accession": [],
        "length": [],
        "e_value": [],
        "score": [],
        "bitscore": [],
        "percent_sim": [],
        "percent_id": [],
    }

    def run(self) -> None:
        """Iterate over MITE fasta files to prepare metadata"""

        for idx, fasta in enumerate(self.fasta.iterdir()):
            with open(fasta) as infile:
                lines = infile.read()
                split_lines = lines.splitlines()
                mite_acc = split_lines[0].split()[0].strip(">")

            self.extract_mite(idx=idx, acc=mite_acc)
            self.extract_xml(acc=mite_acc)

            self.fasta_efi_est.append(f"{lines}\n")

        self.output.mkdir(exist_ok=True)

        df1 = pd.DataFrame(self.metadata_efi_est)
        df1.to_csv(Path(self.output).joinpath("mite_ssn_metadata.csv"), index=False)

        df2 = pd.DataFrame(self.nr_blast_matches)
        df2.to_csv(Path(self.output).joinpath("blast_matches_details.csv"), index=False)

        with open(
            Path(self.output).joinpath("mite_ssn_seqs.fasta"), "w", encoding="utf-8"
        ) as outfile:
            outfile.writelines(self.fasta_efi_est)

    def extract_mite(self, acc: str, idx: int) -> None:
        """Extracts metadata for SSN from mite files

        Arguments:
            acc: a MITE accession
            idx: the current loop index
        """
        with open(self.data.joinpath(f"{acc}.json")) as mite_file:
            mite_data = json.load(mite_file)

        self.metadata_efi_est["key"].append(f"{idx}".rjust(7, "z"))
        self.metadata_efi_est["mite_acc"].append(mite_data["accession"])
        self.metadata_efi_est["enzyme_name"].append(
            mite_data["enzyme"].get("name", "").replace(",", "")
        )
        self.metadata_efi_est["enzyme_description"].append(
            mite_data["enzyme"].get("description", "").replace(",", "")
        )
        self.metadata_efi_est["tailoring"].append(
            "|".join(
                sorted(
                    {
                        tailoring
                        for reaction in mite_data.get("reactions")
                        for tailoring in reaction.get("tailoring", [])
                    }
                )
            )
        )
        self.metadata_efi_est["id_uniprot"].append(
            mite_data["enzyme"]["databaseIds"].get("uniprot", "")
        )
        self.metadata_efi_est["id_genpept"].append(
            mite_data["enzyme"]["databaseIds"].get("genpept", "")
        )
        self.metadata_efi_est["id_mibig"].append(
            mite_data["enzyme"]["databaseIds"].get("mibig", "")
        )

    def extract_xml(self, acc: str) -> None:
        """Extracts metadata for SSN from BLAST XML file

        Counts matches >= 95% similarity and adds to efi-est metadata
        Also stores matches >= 95% similarity for dumping as csv

        Arguments:
            acc: a MITE accession
        """
        with open(self.ncbi_results.joinpath(f"{acc}.xml")) as xml_file:
            blast_record = NCBIXML.read(xml_file)

        counter = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sim_perc = round((hsp.positives / hsp.align_length) * 100, 2)
                id_perc = round((hsp.identities / hsp.align_length) * 100, 2)

                if sim_perc >= 70:
                    self.nr_blast_matches["mite_acc"].append(acc)
                    self.nr_blast_matches["accession"].append(alignment.accession)
                    self.nr_blast_matches["length"].append(alignment.length)
                    self.nr_blast_matches["e_value"].append(hsp.expect)
                    self.nr_blast_matches["score"].append(hsp.score)
                    self.nr_blast_matches["bitscore"].append(hsp.bits)
                    self.nr_blast_matches["percent_sim"].append(sim_perc)
                    self.nr_blast_matches["percent_id"].append(id_perc)

                    counter += 1

        self.metadata_efi_est["ncbi_nr_matches"].append(counter)


class PlotManager(AbstractManager):
    """Organizes code for plotting"""

    def run(self) -> None:
        """Creates boxplot of NCBI NR matches"""

        infile = self.output.joinpath("mite_ssn_metadata.csv")

        if not infile.exists():
            logger.warning(f"Could not find input data '{infile}' - SKIP")

        df = pd.read_csv(infile)
        df.sort_values(by=["ncbi_nr_matches"], ascending=True, inplace=True)

        log_counts = np.log10(df["ncbi_nr_matches"])
        bins = np.linspace(np.floor(min(log_counts)), np.ceil(max(log_counts)), 15)
        hist_vals, bin_edges = np.histogram(log_counts, bins=bins)

        fig, ax = plt.subplots(figsize=(4, 4))
        for i in range(len(hist_vals)):
            ax.bar(bin_edges[i], hist_vals[i],
                   width=bin_edges[i + 1] - bin_edges[i],
                   align='edge',
                   color="darkgrey",
                   linewidth=1,
                   edgecolor='black')

        ax.set_xticks([0, 1, 2, 3, 4])
        ax.set_xticklabels(['1', '10', '100', '1,000', '10,000'])
        ax.set_xlabel('Match counts (log-scale)', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)

        sns.despine()
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.rcParams['font.family'] = 'sans-serif'

        plt.tight_layout()

        plt.savefig(self.output.joinpath("histogram.svg"), format="svg")

        self.write_summary(df)

    def write_summary(self, df: pd.DataFrame) -> None:
        """Write data summary of boxplot"""
        summary = df["ncbi_nr_matches"].describe()

        summary = {
            "count": summary["count"],
            "mean": summary["mean"],
            "std": summary["std"],
            "25%": summary["25%"],
            "50%": summary["50%"],
            "75%": summary["75%"],
            "max": summary["max"],
            "iqr": summary["75%"] - summary["25%"],
            "lower_bound": summary["25%"] - (1.5 * (summary["75%"] - summary["25%"])),
            "upper_bound": summary["75%"] + (1.5 * (summary["75%"] - summary["25%"])),
        }

        json_dict = {key: float(value) for key, value in summary.items()}

        with open(self.output.joinpath("summary_stats.json"), "w") as outfile:
            outfile.write(json.dumps(json_dict, indent=2))


def main() -> None:
    """Function to execute main body of code"""

    parser = argparse.ArgumentParser(description="Script to create MITE SSN files")
    parser.add_argument(
        "-r",
        "--record",
        type=str,
        required=True,
        help="Specifies the Zenodo record ID of the mite_data version to use.",
    )
    args = parser.parse_args(sys.argv[1:])

    download_manager = DownloadManager(record=args.record)
    download_manager.run()

    blast_manager = BlastManager()
    blast_manager.run()

    metadata_manager = MetadataManager()
    metadata_manager.run()

    plot_manager = PlotManager()
    plot_manager.run()


if __name__ == "__main__":
    main()
