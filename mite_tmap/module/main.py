"""Prepare tmap for MITE dataset

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

from annoy import AnnoyIndex
import argparse
from drfp import DrfpEncoder
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from pydantic import BaseModel
import requests
import seaborn as sns
from scipy.spatial.distance import cosine as cosine_distance
import tmap as tm

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)


class AbstractManager(BaseModel):
    """Class to contain general attributes

    data: path to mite_data json files
    output: resulting files for SSN
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
    """Class to generate data file for plotting

    mite_acc: a dict containing mite data
    """

    mite_data: dict = {
        "mite_acc": [],
        "reaction": [],
        "enzyme_name": [],
        "enzyme_description": [],
        "tailoring": [],
        "id_uniprot": [],
        "id_genpept": [],
    }

    def run(self) -> None:
        """Iterates over MITE files to construct input data for plotting"""
        for mite in self.data.iterdir():
            with open(mite) as infile:
                mite_data = json.load(infile)

            if mite_data["status"] != "active":
                logger.info(f"{mite_data['accession']} is retired - SKIP")
                continue

            self.extract_mite(mite_data)

        self.output.mkdir(exist_ok=True)

        df1 = pd.DataFrame(self.mite_data)
        df1.sort_values(by=["tailoring"], ascending=True, inplace=True)
        df1.to_csv(Path(self.output).joinpath("mite_data.csv"), index=False)

    def extract_mite(self, data: dict) -> None:
        """Extract data for tmap plotting

        Onle one representative reaction per MITE entry is plotted

        Arguments:
            data: a dict from mite json file
        """
        categ = "|".join(
            sorted(
                {
                    tailoring
                    for reaction in data.get("reactions")
                    for tailoring in reaction.get("tailoring", [])
                }
            )
        )
        if re.search(r"\|", categ):
            categ = "Multiple"

        self.mite_data["reaction"].append(
            f"{data['reactions'][0]['reactions'][0]['substrate']}>>{'.'.join(data['reactions'][0]['reactions'][0]['products'])}"
        )
        self.mite_data["mite_acc"].append(data["accession"])
        self.mite_data["enzyme_name"].append(
            data["enzyme"].get("name", "").replace(",", "")
        )
        self.mite_data["enzyme_description"].append(
            data["enzyme"].get("description", "").replace(",", "")
        )
        self.mite_data["id_uniprot"].append(
            data["enzyme"]["databaseIds"].get("uniprot", "")
        )
        self.mite_data["id_genpept"].append(
            data["enzyme"]["databaseIds"].get("genpept", "")
        )
        self.mite_data["tailoring"].append(categ)


class PlotManager(AbstractManager):
    """Organizes code to run the tmap plotting"""

    def run(self) -> None:
        """Runs the code"""

        if not self.output.joinpath("rxn_smiles.pickle").exists():
            self.generate_fps()

        self.plot_fps_sns()

    def generate_fps(self) -> None:
        """Generates drfp fingerprints and dumps them as pickle file"""
        df = pd.read_csv(self.output.joinpath("mite_data.csv"))
        rxn_smiles = df["reaction"].tolist()

        enc = DrfpEncoder.encode(rxn_smiles)

        with open(self.output.joinpath("rxn_smiles.pickle"), "wb") as outfile:
            pickle.dump(obj=enc, file=outfile)

    def plot_fps_sns(self) -> None:
        """Reads fingerprints from pickle file and plots them using tmap and seaborn"""
        df = pd.read_csv(self.output.joinpath("mite_data.csv"))

        with open(self.output.joinpath("rxn_smiles.pickle"), "rb") as infile:
            mite_drfp = pickle.load(infile)

        fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)
        fig.set_figheight(5)
        fig.set_figwidth(10)

        CFG_TMAP = tm.LayoutConfiguration()
        CFG_TMAP.k = 50
        CFG_TMAP.kc = 50
        CFG_TMAP.sl_scaling_min = 1.0
        CFG_TMAP.sl_scaling_max = 1.0
        CFG_TMAP.sl_repeats = 1
        CFG_TMAP.sl_extra_scaling_steps = 2
        CFG_TMAP.placer = tm.Placer.Barycenter
        CFG_TMAP.merger = tm.Merger.LocalBiconnected
        CFG_TMAP.merger_factor = 2.0
        CFG_TMAP.merger_adjustment = 0
        CFG_TMAP.fme_iterations = 1000
        CFG_TMAP.sl_scaling_type = tm.ScalingType.RelativeToDesiredLength
        CFG_TMAP.node_size = 5
        CFG_TMAP.mmm_repeats = 1

        y_values = df["tailoring"].tolist()
        knn = []

        annoy = AnnoyIndex(2048, metric="angular")
        for i, v in enumerate(mite_drfp):
            annoy.add_item(i, v)

        annoy.build(10)
        for i in range(len(mite_drfp)):
            for j in annoy.get_nns_by_item(i, 10):
                knn.append((i, j, cosine_distance(mite_drfp[i], mite_drfp[j])))

        x, y, s, t, _ = tm.layout_from_edge_list(len(mite_drfp), knn, config=CFG_TMAP)
        for i in range(len(s)):
            ax.plot(
                [x[s[i]], x[t[i]]],
                [y[s[i]], y[t[i]]],
                "k-",
                linewidth=0.5,
                alpha=0.5,
                zorder=1,
            )

        df_tmap = pd.DataFrame({"x": x, "y": y, "c": y_values})

        # full palette
        # sns.set_palette(
        #     [
        #         "#9ad59a",  # acetylation
        #         "#91efd4",  # amination
        #         "#c2c3d0",  # biaryl
        #         "#ccc0dd",  # cycliz
        #         "#decd87",  # deamination
        #         "#dcb89d",  # decarboxylation
        #         "#ff9b9b",  # dehydr
        #         "#ffaf79",  # dehydrogen
        #         "#ffffaf",  # dioxygenation
        #         "#78804d",  # epimerisation
        #         "#c0d1b6",  # glycosly
        #         "#93b0b0",  # halogenation
        #         "#8caad1",  # heterocycl
        #         "#9d6ab1",  # hydrolysis
        #         "#b958aa",  # hydoxylation
        #         "#ab4962",  # macrolactam
        #         "#dd626a",  # methylation
        #         "#d57058",  # monooxygenation
        #         "#5599ff",  # multiple
        #         "#afafaf",  # other
        #         "#b38466",  # oxidation
        #         "#5dc851",  # prenylation
        #         "#958479",  # reduction
        #         "#f78c6a",  # sulfonation
        #     ]
        # )
        # simplified palette
        sns.set_palette(
            [
                "#afafaf",  # acetylation
                "#afafaf",  # amination
                "#afafaf",  # biaryl
                "#afafaf",  # cycliz
                "#afafaf",  # deamination
                "#afafaf",  # decarboxylation
                "#afafaf",  # dehydr
                "#afafaf",  # dehydrogen
                "#afafaf",  # dioxygenation
                "#afafaf",  # epimerisation
                "#afafaf",  # glycosly
                "#5599ff",  # halogenation
                "#afafaf",  # heterocycl
                "#afafaf",  # hydrolysis
                "#b958aa",  # hydoxylation
                "#afafaf",  # macrolactam
                "#dd626a",  # methylation
                "#afafaf",  # monooxygenation
                "#afafaf",  # multiple
                "#afafaf",  # other
                "#b38466",  # oxidation
                "#afafaf",  # prenylation
                "#afafaf",  # reduction
                "#afafaf",  # sulfonation
            ]
        )

        sns.scatterplot(
            x="x",
            y="y",
            hue="c",
            data=df_tmap,
            palette=sns.color_palette(),
            s=75.0,
            ax=ax,
            zorder=2,
        )

        legend = ax.legend(
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            fancybox=False,
            shadow=False,
            frameon=False,
            ncol=1,
            fontsize=8,
        )

        ax.axis("off")
        plt.tight_layout()

        plt.savefig(self.output.joinpath("tmap_mite_plain.svg"), format="svg")

        for i in range(len(df_tmap)):
            plt.text(
                x=df_tmap["x"][i],
                y=df_tmap["y"][i],
                s=str(df["enzyme_name"][i]),
                fontsize=6,
                ha="center",
            )

        plt.savefig(self.output.joinpath("tmap_mite_annotated.svg"), format="svg")


def main() -> None:
    """Function to execute main body of code"""

    parser = argparse.ArgumentParser(
        description="Script to create tmap plot for MITE dataset"
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

    plot_manager = PlotManager()
    plot_manager.run()


if __name__ == "__main__":
    main()
