# Copyright (c) cifcomb Development Team.
# Distributed under the terms of the MIT License.

import os
import numpy as np
import itertools

import argparse
import warnings

from pymatgen.core import Composition

import cifcombo.core as core
import cifcombo.const as const

warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-md", "--makedatabase", help="-m cifdir1 cifdir2 ...", nargs="*"
    )
    parser.add_argument(
        "-wc",
        "--withoutcifcopy",
        help="make a database without copying the cif dirs (just make symbolic links to them).",
        action="store_true",
    )
    parser.add_argument(
        "-lc",
        "--makeciflink",
        help="-lc cifdir, make a link to the cif dir (without copy)",
        type=str,
    )
    parser.add_argument(
        "-ld",
        "--makedatalink",
        help="-ld datadir, make a link to the data dir (without copy)",
        type=str,
    )
    parser.add_argument(
        "-cc", "--copycifdir", help="-cc cifdir, copy the cif dir", type=str
    )
    parser.add_argument(
        "-cd",
        "--copydatadir",
        help="-cd datadir, copy the data dir",
        type=str,
    )
    parser.add_argument(
        "-s",
        "--searchcombination",
        help="-s composition1 composition2 ... (Search for combinations of the given composition from the database. When you specify -s, cifcombo tries to decompose the target composition with the given compositions.)",
        nargs="*",
    )
    parser.add_argument(
        "-t",
        "--targetcompsition",
        help="-t targetcomposition (When you specify -s, cifcombo tries to decompose the target composition with the given compositions. When you do not specify -s, cifcombo searches possible decompositions (the default num. = 2, you can specify the number by -nd) from the database.)",
        type=str,
    )
    parser.add_argument(
        "-nd",
        "--numpossibledecomposition",
        help="-nd the num. of possible decompositions valid only with -t and without -s",
        type=int,
        default=2,
    )
    parser.add_argument("-g", "--getcif", help="-g [cifid] ...", nargs="*")

    args = parser.parse_args()

    if args.makedatabase:
        core.make_database(
            cif_dir_list=args.makedatabase, cif_dir_link=args.withoutcifcopy
        )

    if args.makeciflink:
        src_cif_dir = const.CIF_DIR
        dst_cif_dir = os.path.abspath(args.makeciflink)
        core.make_link(
            src=os.path.relpath(src_cif_dir), dst=os.path.relpath(dst_cif_dir)
        )

    if args.makedatalink:
        src_data_dir = const.DATA_DIR
        dst_data_dir = os.path.abspath(args.makedatalink)
        core.make_link(
            src=os.path.relpath(src_data_dir),
            dst=os.path.relpath(dst_data_dir),
        )

    if args.copycifdir:
        src_cif_dir = os.path.abspath(args.copycifdir)
        dst_cif_dir = const.CIF_DIR
        core.copy_dir(
            src=os.path.relpath(src_cif_dir), dst=os.path.relpath(dst_cif_dir)
        )

    if args.copydatadir:
        src_data_dir = os.path.abspath(args.copydatadir)
        dst_data_dir = const.DATA_DIR
        core.copy_dir(
            src=os.path.relpath(src_data_dir),
            dst=os.path.relpath(dst_data_dir),
        )

    if args.searchcombination:
        compound_list = args.searchcombination

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

                if len(input_nominal_ratio_list) > 1:
                    print(
                        f"{target_nominal_ratio:.0f}*{args.targetcompsition} = {eq}"
                    )

                else:
                    print(f"{eq}")

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

                    if len(input_compositions_list) > 1:
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
                        eq = " + ".join(
                            [
                                f"{str(input_compositions_list[i]).replace(' ', '')}"
                                for i in range(len(input_compositions_list))
                            ]
                        )
                        print(f"{eq}, cifid={cifid}, SG={space_group}")

            else:
                print("No cif is found in the cif database.")
                print("If you try to synthesize a new compound,")
                print("plz. specify the composition with -t option.")

    else:
        if args.targetcompsition:
            target_composition = Composition(args.targetcompsition)
            print(f"target composition = {target_composition}")
            candidate_el_amt_dict = (
                target_composition.element_composition.get_el_amt_dict()
            )
            unique_element_list = candidate_el_amt_dict.keys()
            # print(unique_element_list)
            all_possible_combinations = []
            for num_c in range(len(unique_element_list)):
                all_possible_combinations += list(
                    itertools.combinations(unique_element_list, num_c + 1)
                )

            # print(all_possible_combinations)
            print(
                f"The num. possible decomposition = {args.numpossibledecomposition}"
            )
            all_possible_compositions = (
                itertools.combinations_with_replacement(
                    all_possible_combinations, args.numpossibledecomposition
                )
            )
            print(
                "Searching for possible decompositions of the target composition..."
            )

            print(
                "It may take a copuple of minites... Plz. be patient ... :-)"
            )

            # print(list(all_possible_compositions))

            for possible_compositions in all_possible_compositions:
                element_included_list = list(
                    set(itertools.chain.from_iterable(possible_compositions))
                )

                if all(
                    [
                        element in element_included_list
                        for element in unique_element_list
                    ]
                ):
                    # print(f"possible_composition={possible_compositions}")
                    all_product_ingredients = []
                    for compositions in possible_compositions:
                        # print(f"compositions = {list(compositions)}")
                        (
                            return_candidate_nominal_ratio_list,
                            return_candidate_composition_list,
                            return_input_nominal_ratio_list,
                            return_input_compositions_list,
                            return_cifid_list,
                            return_space_group_list,
                        ) = core.search_combination(
                            compounds_list=compositions
                        )

                        # print(set(return_candidate_composition_list))

                        all_product_ingredients.append(
                            set(return_candidate_composition_list)
                        )

                    p = itertools.product(*all_product_ingredients)

                    for input_composition_list in p:
                        # print(input_composition_list)
                        (
                            target_nominal_ratio,
                            input_nominal_ratio_list,
                        ) = core.decompose_composition(
                            target_composition=target_composition,
                            input_composition_list=input_composition_list,
                        )

                        if target_nominal_ratio is not np.nan:
                            eq = " + ".join(
                                [
                                    f"{input_nominal_ratio_list[i]:.0f}*{str(input_composition_list[i]).replace(' ', '')}"
                                    for i in range(len(input_composition_list))
                                ]
                            )

                            print(
                                f"  {target_nominal_ratio:.0f}*{args.targetcompsition} = {eq}"
                            )
            print("Done.")

    if args.getcif:
        cifid_list = args.getcif
        core.getcif(cifid_list=cifid_list)


if __name__ == "__main__":
    main()
