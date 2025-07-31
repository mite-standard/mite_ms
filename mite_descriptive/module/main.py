"""Prepare descriptive plots for MITE dataset

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
import re
from pathlib import Path
import shutil
import sys

import argparse
import matplotlib.pyplot as plt
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
    output: resulting files for SSN
    """

    data: Path = Path(__file__).parent.joinpath("data/data")
    fasta: Path = Path(__file__).parent.joinpath("data/fasta")
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
                "Could not determine data storage location in downloaded directory."
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


class MetadataManager(AbstractManager):
    """Class to generate metadata file

    metadata: a dict containing metadata
    datapoints: all mite entries flattened to count datapoints
    """

    metadata: dict = {
        "key": [],
        "mite_acc": [],
        "enzyme_name": [],
        "enzyme_description": [],
        "tailoring": [],
        "id_uniprot": [],
        "id_genpept": [],
        "id_mibig": [],
        "rhea_mibig": [],
    }
    datapoints: list = []

    def run(self) -> None:
        """Iterate over MITE fasta files to prepare metadata"""

        for idx, fasta in enumerate(self.fasta.iterdir()):
            with open(fasta) as infile:
                lines = infile.read()
                split_lines = lines.splitlines()
                mite_acc = split_lines[0].split()[0].strip(">")

            self.extract_mite(idx=idx, acc=mite_acc)

        self.output.mkdir(exist_ok=True)

        df1 = pd.DataFrame(self.metadata)
        df1.to_csv(Path(self.output).joinpath("mite_metadata.csv"), index=False)

        with open(Path(self.output).joinpath("flat_mite.csv"), 'w') as f:

            cleaned = []
            for line in self.datapoints:
                if line[0].startswith("changelog"):
                    continue
                else:
                    cleaned.append(line[1])

            for line in set(cleaned):
                f.write(f"{line}\n")

    def extract_mite(self, acc: str, idx: int) -> None:
        """Extracts metadata for SSN from mite files

        Arguments:
            acc: a MITE accession
            idx: the current loop index
        """
        with open(self.data.joinpath(f"{acc}.json")) as mite_file:
            mite_data = json.load(mite_file)

        flattened = self.flatten(mite_data)
        self.datapoints.extend(flattened)

        self.metadata["key"].append(f"{idx}".rjust(7, "z"))
        self.metadata["mite_acc"].append(mite_data["accession"])
        self.metadata["enzyme_name"].append(
            mite_data["enzyme"].get("name", "").replace(",", "")
        )
        self.metadata["enzyme_description"].append(
            mite_data["enzyme"].get("description", "").replace(",", "")
        )
        categ = "|".join(
            sorted(
                {
                    tailoring
                    for reaction in mite_data.get("reactions")
                    for tailoring in reaction.get("tailoring", [])
                }
            )
        )
        if re.search(r"\|", categ):
            self.metadata["tailoring"].append("Multiple")
        else:
            self.metadata["tailoring"].append(categ)

        self.metadata["id_uniprot"].append(
            mite_data["enzyme"]["databaseIds"].get("uniprot", "")
        )
        self.metadata["id_genpept"].append(
            mite_data["enzyme"]["databaseIds"].get("genpept", "")
        )
        self.metadata["id_mibig"].append(
            mite_data["enzyme"]["databaseIds"].get("mibig", "")
        )

        flag_rhea = False
        try:
            for reaction in mite_data["reactions"]:
                if reaction.get("databaseIds", {}).get("rhea"):
                    flag_rhea = True
        except KeyError:
            logger.info(f"{acc} does not have a database reference - pass")

        if flag_rhea and mite_data["enzyme"]["databaseIds"].get("mibig"):
            self.metadata["rhea_mibig"].append("Rhea+MIBiG")
        elif flag_rhea:
            self.metadata["rhea_mibig"].append("Rhea")
        elif mite_data["enzyme"]["databaseIds"].get("mibig"):
            self.metadata["rhea_mibig"].append("MIBiG")
        else:
            self.metadata["rhea_mibig"].append("No crosslink")

    def flatten(self, data: dict, parent_key=''):
        """Extracts data points from MITE entries

        Args:
            data: a MITE dict
        """
        items = []

        if isinstance(data, dict):
            for k, v in data.items():
                new_key = f"{parent_key}.{k}" if parent_key else k
                items.extend(self.flatten(v, new_key))
        elif isinstance(data, list):
            for i, v in enumerate(data):
                new_key = f"{parent_key}[{i}]"
                items.extend(self.flatten(v, new_key))
        else:
            items.append((parent_key, data))

        return items



class PlotManager(AbstractManager):
    """Organizes code for plotting"""

    def run(self) -> None:

        infile = self.output.joinpath("mite_metadata.csv")

        if not infile.exists():
            logger.warning(f"Could not find input data '{infile}' - SKIP")

        df = pd.read_csv(infile)

        self.plot_countplot(df)
        self.plot_pieplot(df)


    def plot_countplot(self, df: pd.DataFrame):
        """Plot of tailoring labels"""

        df_sorted = df.sort_values(by=["tailoring"], ascending=True)
        plt.figure(figsize=(8, 6))
        sns.countplot(y="tailoring", data=df_sorted)

        plt.tight_layout()
        plt.xlabel("Category", fontsize=10)

        plt.savefig(self.output.joinpath("countplot.svg"), format="svg")
        plt.clf()
        plt.close()

    def plot_pieplot(self, df: pd.DataFrame):
        """Plot of mibig/rhea crosslink labels"""
        df_sorted = df.sort_values(by=["rhea_mibig"], ascending=True)
        value_counts = df_sorted["rhea_mibig"].value_counts()

        plt.figure(figsize=(4, 4))
        plt.pie(
            value_counts,
            labels=value_counts.index,
            autopct="%1.1f%%",
            startangle=90,
            colors=sns.color_palette("Set2", len(value_counts)),
        )

        plt.tight_layout()

        plt.savefig(self.output.joinpath("pieplot.svg"), format="svg")
        plt.clf()
        plt.close()


def main() -> None:
    """Function to execute main body of code"""

    parser = argparse.ArgumentParser(
        description="Script to create descriptive plots for MITE dataset"
    )
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

    metadata_manager = MetadataManager()
    metadata_manager.run()

    plot_manager = PlotManager()
    plot_manager.run()


if __name__ == "__main__":
    main()
