
from pathlib import Path

import argparse
import pandas as pd
import plotly.express as px
import requests

src = Path(__file__).parent.joinpath("summary.csv")


def download() -> pd.DataFrame:
    """Downloads source file from https://mite.bioinformatics.nl/downloads/mite_overview"""

    if src.exists():
        return

    response = requests.get(
        "https://mite.bioinformatics.nl/downloads/mite_overview"
    )
    if response.status_code != 200:
        raise RuntimeError("Error fetching data")

    with open(src, "w") as f:
        f.write(response.text)

def plot() -> None:
    """Plot sunburst from data file"""

    inner = "kingdom"
    outer = "phylum"

    df = pd.read_csv(src)
    df_clean = df[(df[inner] != "Not found") & (df[outer] != "Not found")]
    sunburst_data = df_clean.groupby([inner, outer]).size().reset_index(name='count')

    print("Number of retained entries for sunburst plot:")
    print(sunburst_data["count"].sum())

    fig = px.sunburst(
        sunburst_data,
        path=[inner, outer],
        values="count",
        color=inner,
    )

    fig.update_layout(
        width=250,
        height=250,
        margin=dict(t=30, l=30, r=30, b=30),
        font=dict(size=12),
        uniformtext=dict(minsize=12, mode='show')
    )

    fig.write_image("mite_sunburst_plot.svg", width=250, height=250)


def main() -> None:
    parser = argparse.ArgumentParser(description="Script to create MITE sunburst plot")
    parser.parse_args()
    download()
    plot()


if __name__ == "__main__":
    main()
