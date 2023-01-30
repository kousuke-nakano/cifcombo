# Copyright (c) cifcomb Development Team.
# Distributed under the terms of the MIT License.

import os
import numpy as np

import argparse
import warnings

from pymatgen.core import Composition

import cifcombo.core as core

warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-m", "--makedatabase", help="-m cifdir1 cifdir2 ...", nargs="*"
    )
    parser.add_argument(
        "-link",
        "--makeciflink",
        help="make links for the cif dirs (without copy)",
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--searchdatabase",
        help="-s composition1 composition2 ...",
        nargs="*",
    )
    parser.add_argument(
        "-t", "--targetcompsition", help="-t targetcomposition ...", type=str
    )
    parser.add_argument("-g", "--getcif", help="-g [cifid] ...", nargs="*")

    args = parser.parse_args()

    if args.makedatabase:
        core.make_database(
            cif_dir_list=args.makedatabase, cif_dir_link=args.makeciflink
        )

    if args.searchdatabase:
        compound_list = args.searchdatabase

        if args.targetcompsition:
            (
                target_nominal_ratio,
                input_nominal_ratio_list,
            ) = core.decompose_composition(
                target_composition=Composition(args.targetcompsition),
                input_composition_list=[Composition(c) for c in compound_list],
            )

            if target_nominal_ratio is not np.nan:
                eq = " + ".join(
                    [
                        f"{input_nominal_ratio_list[i]:.0f}*{str(compound_list[i]).replace(' ', '')}"
                        for i in range(len(compound_list))
                    ]
                )

                print(
                    f"{target_nominal_ratio:.0f}*{args.targetcompsition} = {eq}"
                )
            else:
                print("No solution is found.")

        else:
            (
                return_candidate_nominal_ratio_list,
                return_candidate_composition_list,
                return_input_nominal_ratio_list,
                return_input_compositions_list,
                return_cifid_list,
                return_space_group_list,
            ) = core.search_combination(compounds_list=compound_list)

            if len(return_candidate_nominal_ratio_list) != 0:
                for (
                    candidate_nominal_ratio,
                    candidate_composition,
                    input_nominal_ratio_list,
                    input_compositions_list,
                    cifid,
                    space_group,
                ) in zip(
                    return_candidate_nominal_ratio_list,
                    return_candidate_composition_list,
                    return_input_nominal_ratio_list,
                    return_input_compositions_list,
                    return_cifid_list,
                    return_space_group_list,
                ):

                    eq = " + ".join(
                        [
                            f"{input_nominal_ratio_list[i]:.0f}*{str(input_compositions_list[i]).replace(' ', '')}"
                            for i in range(len(input_compositions_list))
                        ]
                    )

                    print(
                        f"{candidate_nominal_ratio:.0f}*{candidate_composition} = {eq}, cifid={cifid}, SG={space_group}"
                    )
            else:
                print("No cif is found in the cif database.")
                print("If you try to synthesize a new compound,")
                print("plz. specify the composition with -t option.")

    if args.getcif:
        cifid_list = args.getcif
        core.getcif(cifid_list=cifid_list)


if __name__ == "__main__":
    main()
