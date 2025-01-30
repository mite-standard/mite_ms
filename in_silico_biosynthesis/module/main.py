"""Proof-of-concept using mite_data for in silico biosynthesis: bottromycin A2 biosynthetic pathway

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


from rdkit import RDLogger
from rdkit.Chem import CanonSmiles, MolFromSequence, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
import requests

RDLogger.DisableLog("rdApp.*")

def apply_reaction(substrate: str, reaction: str) -> str:
    """Applies reaction SMARTS to substrate and returns the product smiles

    Args:
        substrate: a SMILES string
        reaction: a reaction SMARTS string

    Returns:
        The product SMILES string
    """
    rd_substrate = MolFromSmiles(substrate)
    rd_reaction = ReactionFromSmarts(reaction)
    products = rd_reaction.RunReactants([rd_substrate])
    products = (MolToSmiles(product[0]) for product in products)

    return next(iter(products))

def mite_api_call(accession: str) -> str:
    """Retrieve json data from MITE API

    Attributes:
        accession: a MITE accession ID

    Raises:
        RuntimeError: Could not connect successfully to MITE API

    Returns:
        the reaction SMARTS of the first reaction entry
    """
    response = requests.get(f"https://mite.bioinformatics.nl/api/{accession}", timeout=10)
    if response.status_code != 200:
        raise RuntimeError(f"Error in MITE API call: {response.status_code}")

    mite_data = response.json()
    return mite_data["reactions"][0]["reactionSMARTS"]


def main() -> None:
    """Run the bottromycin in silico biosynthesis"""

    synthesis = {
        "MITE0000048": {"enzyme": "BotP"},
        "MITE0000049": {"enzyme": "botRMT1"},
        "MITE0000046": {"enzyme": "botRMT2"},
        "MITE0000047": {"enzyme": "botRMT3"},
        "MITE0000052": {"enzyme": "BotC"},
        "MITE0000053": {"enzyme": "BotCD"},
        "MITE0000050": {"enzyme": "BotAH"},
        "MITE0000051": {"enzyme": "BotH"},
        "MITE0000044": {"enzyme": "BotCYP"},
        "MITE0000045": {"enzyme": "BotOMT"},
    }

    for key in synthesis.keys():
        synthesis[key]["r_smarts"] = mite_api_call(accession=key)

    products = [MolToSmiles(MolFromSequence("MGPVVVFDCMTADFLNDDPNNAELSALEMEELESWGAWDGEATS"))] # bottromycin A2 precursor from Streptomyces sp. BC16019 (MIBiG BGC0000469)

    for key, value in synthesis.items():
        print(f"### Applying reaction {key}")
        try:
            product = apply_reaction(substrate=products[-1], reaction=value["r_smarts"])
            products.append(product)
            print(f"###### Resulting product: {product}")
        except StopIteration:
            print("### The last reaction SMARTS failed to yield any product.")
            print(f"### The most recent product is: {products[-1]}")
            exit(1)

    gen_bottromycin_a2 = CanonSmiles(products[-1])
    real_bottromycin_a2 = CanonSmiles(
        "C[C@H]1[C@H]2C(N[C@@H](C(C)C)C(N[C@@H](C(C)(C)C)/C(=N/[C@@H](C(C)(C)C)C(N[C@H](C(N[C@@H](c3sccn3)CC(OC)=O)=O)[C@H](c3ccccc3)C)=O)/NCC(=O)N2CC1)=O)=O"
    )  # as described by doi:10.1039/D0NP00097C

    print(f"In silico generated bottromycin A2: {gen_bottromycin_a2}")
    print(f"Literature-reported bottromycin A2: {real_bottromycin_a2}")

    if gen_bottromycin_a2 == real_bottromycin_a2:
        print(
            "SUCCESS: In silico generated bottromycin A2 corresponds to literature-reported one."
        )
    else:
        print("FAILURE: In silico generated bottromycin A2 DOES NOT correspond to literature-reported one")

if __name__ == "__main__":
    main()