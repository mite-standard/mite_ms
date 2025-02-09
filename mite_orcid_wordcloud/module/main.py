"""Prepare wordcloud plot for MITE dataset

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

import argparse
import pandas as pd
from pydantic import BaseModel
import requests


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

class AbstractManager(BaseModel):
    """Class to contain general attributes

    data: path to mite_data json files
    output: output folder
    """

    data: Path = Path(__file__).parent.joinpath("data/data")
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

        os.remove(self.record)
        shutil.rmtree(self.record_unzip)

class DataManager(AbstractManager):
    """Class to collect data and plot wordcloud

    orcids: a list of orcids
    """

    orcids: list = []

    def run(self) -> None:
        """Iterate over MITE fasta files to prepare metadata"""

        self.extract_data()

    def extract_data(self) -> None:
        """Pull out orcids from mite files and dump it"""

        for mite in self.data.iterdir():
            with open(mite) as infile:
                mite_data = json.load(infile)

            for log in mite_data["changelog"]:
                self.orcids.extend(log["contributors"])
                self.orcids.extend(log["reviewers"])

        self.orcids = [word for word in self.orcids if not word == "AAAAAAAAAAAAAAAAAAAAAAAA"]

        self.output.mkdir(exist_ok=True)

        frequency = {}
        for word in self.orcids:
            if word not in frequency:
                frequency[word] = 1
            else:
                frequency[word] += 1

        frame = {
            "weight": [value for value in frequency.values()],
            "word": [key for key in frequency.keys()]
        }

        df = pd.DataFrame(frame)
        df.to_csv(self.output.joinpath("orcids.csv"), index=False)

def main() -> None:
    """Function to execute main body of code"""

    parser = argparse.ArgumentParser(
        description="Script to create ORCID wordcloud for MITE dataset"
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

    data_manager = DataManager()
    data_manager.run()


if __name__ == "__main__":
    main()

