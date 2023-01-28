# Copyright (c) cifcomb Development Team.
# Distributed under the terms of the MIT License.

import os


import argparse
import warnings

import core

warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"


if __name__ == "__main__":

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
        "-s", "--searchdatabase", help="-s compound1 compound2 ...", nargs="*"
    )
    parser.add_argument("-g", "--getcif", help="-g [cifid] ...", nargs="*")

    args = parser.parse_args()

    if args.makedatabase:
        core.make_database(
            cif_dir_list=args.makedatabase, cif_dir_link=args.makeciflink
        )

    if args.searchdatabase:
        compound_list = args.searchdatabase

        (
            return_candidate_composition_list,
            return_input_nominal_ratio_list,
            return_input_compositions_list,
            return_codid_list,
            return_space_group_list,
        ) = core.search_combination(compounds_list=compound_list)

        for (
            candidate_composition,
            input_nominal_ratio_list,
            input_compositions_list,
            codid,
            space_group,
        ) in zip(
            return_candidate_composition_list,
            return_input_nominal_ratio_list,
            return_input_compositions_list,
            return_codid_list,
            return_space_group_list,
        ):

            eq = " + ".join(
                [
                    f"{input_nominal_ratio_list[i]:.0f}*{str(input_compositions_list[i]).replace(' ', '')}"
                    for i in range(len(input_compositions_list))
                ]
            )

            print(
                f"{candidate_composition} = {eq}, codid={codid}, SG={space_group}"
            )

    if args.getcif:
        cifid_list = args.getcif
        core.getcif(cifid_list=cifid_list)
