# Copyright (c) cifcomb Development Team.
# Distributed under the terms of the MIT License.

import os
import io
import itertools
import numpy as np
from tabulate import tabulate
from fnmatch import fnmatch
from typing import Tuple
import hashlib
import warnings
import joblib
from glob import glob
import pickle
import pandas as pd
import shutil

import logging

from tqdm.auto import tqdm
import contextlib
from typing import Optional

from monty.fractions import gcd_float
from pymatgen.core import Composition, PeriodicSite, Structure
from pymatgen.core.composition import reduce_formula
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import in_coord_list_pbc
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.util.string import formula_double_format

import cifcomb.const as const

warnings.simplefilter("ignore")
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"

INT_THRESHOLD = 1.0e-5
SYM_TOLERANCE = 0.001
ANGLE_TOLERANCE = 5.0


@contextlib.contextmanager
def tqdm_joblib(
    total: Optional[int] = None, miniters: Optional[int] = 1, **kwargs
):

    pbar = tqdm(total=total, miniters=miniters, smoothing=0, **kwargs)

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            pbar.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback

    try:
        yield pbar
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        pbar.close()


def make_database(cif_dir_list=[], cif_dir_link=True):
    os.makedirs(const.DATA_DIR, exist_ok=True)
    os.makedirs(const.CIF_DIR, exist_ok=True)
    if cif_dir_link:
        for cif_dir_path in cif_dir_list:
            cif_path = os.path.join(
                const.CIF_DIR, os.path.basename(cif_dir_path.rstrip("/"))
            )
            if os.path.isfile(cif_path):
                if os.path.islink(cif_path):
                    os.unlink(cif_path)
                else:
                    os.remove(cif_path)
            if os.path.isdir(cif_path):
                if os.path.islink(cif_path):
                    os.unlink(cif_path)
                else:
                    shutil.rmtree(cif_path)
            os.symlink(
                os.path.abspath(cif_dir_path),
                os.path.join(
                    const.CIF_DIR, os.path.basename(cif_dir_path.rstrip("/"))
                ),
                target_is_directory=True,
            )
    else:
        for cif_dir_path in cif_dir_list:
            shutil.copytree(
                cif_dir_path,
                os.path.join(
                    const.CIF_DIR, os.path.basename(cif_dir_path.rstrip("/"))
                ),
                dirs_exist_ok=True,
            )
    cif_list = glob(os.path.join(const.CIF_DIR, "**", "*.cif"), recursive=True)

    col = [
        "cifid",
        "lattice_type",
        "point_group",
        "space_group_num",
        "space_group",
        "chem_formula",
        "sub_chem_formula",
    ]

    num_thread = int(os.cpu_count() / 2)  # -1 -> max available
    print(f"used num_threads = {num_thread}")

    with tqdm_joblib(
        len(cif_list),
        miniters=int(len(cif_list) / const.JOBLIB_MINITERS_DENOM),
        maxinterval=const.JOBLIB_MAX_INTERVALS,
    ):
        results = joblib.Parallel(n_jobs=num_thread)(
            [joblib.delayed(return_chemical_system)(cif) for cif in cif_list]
        )

    record_all_list = []
    record_dict = {}

    for element_key, record in results:
        if element_key == "NA":
            continue
        record_all_list.append(record)
        if element_key in record_dict.keys():
            record_dict[element_key].append(record)
        else:
            record_dict[element_key] = [record]

    record_dict_pickled = {
        hash: pd.DataFrame(record, columns=col)
        for hash, record in record_dict.items()
    }

    with open(os.path.join(const.DATA_DIR, "cif_summary.pkl"), "wb") as f:
        pickle.dump(record_dict_pickled, f)

    pd.DataFrame(record_all_list, columns=col).to_csv(
        os.path.join(const.DATA_DIR, "cif_summary.csv"), index=False
    )

    cif_path_data = {
        os.path.splitext(os.path.basename(cif_path))[0]: os.path.relpath(
            cif_path, const.CIF_DIR
        )
        for cif_path in cif_list
    }
    with open(
        os.path.join(const.DATA_DIR, "cifpath_data_summary.pkl"), "wb"
    ) as f:
        pickle.dump(cif_path_data, f)


def getcif(cifid_list=[]):
    current_dir = os.getcwd()
    with open(
        os.path.join(const.DATA_DIR, "cifpath_data_summary.pkl"), "rb"
    ) as f:
        cif_data_pickled = pickle.load(f)

    for cifid in cifid_list:
        try:
            print(f"copying {cifid}.cif to the current directory.")
            path_to_cif = cif_data_pickled[cifid]
            shutil.copyfile(
                os.path.join(const.CIF_DIR, path_to_cif),
                os.path.join(current_dir, f"{cifid}.cif"),
            )

        except KeyError:
            print("{cifid}.cif is not found")


def decompose_composition(
    target_composition, input_composition_list=[], max_multiplication=10
):
    candidate_el_amt_dict = (
        target_composition.element_composition.get_el_amt_dict()
    )
    unique_element_list = candidate_el_amt_dict.keys()

    """
    n*TiO2 + m*BaO = 42Ba + 32Ti + 68O
    ->
    0n + 1m = 42 (Ba) * k
    1n + 0m = 32 (Ti) * k
    2n + 1m = 68 (O) * k
    where k is an interger
    """

    for iii in range(max_multiplication):
        multi = iii + 1
        A = np.zeros(
            shape=(len(unique_element_list), len(input_composition_list))
        )
        B = np.zeros(shape=(len(unique_element_list), 1))

        for n_element, element in enumerate(unique_element_list):
            B[n_element, 0] = int(candidate_el_amt_dict[element]) * multi

            for n_composition, composition in enumerate(
                input_composition_list
            ):
                input_el_amt_dict = (
                    composition.element_composition.get_el_amt_dict()
                )

                if element in input_el_amt_dict.keys():
                    A[n_element, n_composition] = float(
                        input_el_amt_dict[element]
                    )
                else:
                    A[n_element, n_composition] = float(0)

        # check if the obtained solution is exact or just the least-squared one.
        X = np.dot(np.linalg.pinv(A), B)  # pseudo-inverse matrix
        B_ = np.dot(A, X)
        diff_B = B - B_
        check_interger = [np.abs(b) < INT_THRESHOLD for b in diff_B.flatten()]
        if not all(check_interger):
            # print("Not an exact solution")
            continue

        input_nominal_ratio_list = X.flatten()
        target_nominal_ratio = multi

        # check if all the nominal composition is integer
        check_interger = [
            np.abs(ratio - round(ratio)) < const.INT_THRESHOLD
            and round(ratio) > 0
            for ratio in input_nominal_ratio_list
        ]
        if all(check_interger):
            target_nominal_ratio = round(target_nominal_ratio)
            input_nominal_ratio_list = [
                round(k) for k in input_nominal_ratio_list
            ]
            return target_nominal_ratio, input_nominal_ratio_list

    # if no solution is found. return np.nan.
    target_nominal_ratio = np.nan
    input_nominal_ratio_list = [np.nan for _ in input_composition_list]

    return target_nominal_ratio, input_nominal_ratio_list


def search_combination(compounds_list=[]):
    # e.g., compounds_list = ["CaO", "SiO2", "O2"]
    input_composition_list = [
        Composition(compound) for compound in compounds_list
    ]
    sum_composition = Composition("".join(compounds_list))
    element_key = sum_composition.chemical_system

    with open(os.path.join(const.DATA_DIR, "cif_summary.pkl"), "rb") as f:
        record_dict_pickled = pickle.load(f)

    try:
        data_pd = record_dict_pickled[element_key]
        cifid_list = list(data_pd["cifid"])
        space_group_list = list(data_pd["space_group"])
        sub_chem_formula_list = list(data_pd["sub_chem_formula"])

        return_candidate_nominal_ratio_list = []
        return_candidate_composition_list = []
        return_input_nominal_ratio_list = []
        return_input_compositions_list = []
        return_cifid_list = []
        return_space_group_list = []

        for cifid, space_group, sub_chem_formula in zip(
            cifid_list, space_group_list, sub_chem_formula_list
        ):
            candidate_composition = Composition(sub_chem_formula)

            (
                target_nominal_ratio,
                input_nominal_ratio_list,
            ) = decompose_composition(
                target_composition=candidate_composition,
                input_composition_list=input_composition_list,
            )

            if target_nominal_ratio is not np.nan:
                return_candidate_nominal_ratio_list.append(
                    target_nominal_ratio
                )
                return_candidate_composition_list.append(candidate_composition)
                return_input_nominal_ratio_list.append(
                    input_nominal_ratio_list
                )
                return_input_compositions_list.append(input_composition_list)
                return_cifid_list.append(cifid)
                return_space_group_list.append(space_group)

    except KeyError:
        return_candidate_nominal_ratio_list = []
        return_candidate_composition_list = []
        return_input_nominal_ratio_list = []
        return_input_compositions_list = []
        return_cifid_list = []
        return_space_group_list = []

    finally:
        return (
            return_candidate_nominal_ratio_list,
            return_candidate_composition_list,
            return_input_nominal_ratio_list,
            return_input_compositions_list,
            return_cifid_list,
            return_space_group_list,
        )


def joblib_search_combination(combination_list=[], num_thread=-1):

    with tqdm_joblib(
        len(combination_list),
        miniters=int(len(combination_list) / const.JOBLIB_MINITERS_DENOM),
        maxinterval=const.JOBLIB_MAX_INTERVALS,
        disable=const.DISABLE_PROGRESSBAR,
    ):
        results = joblib.Parallel(n_jobs=num_thread)(
            [
                joblib.delayed(search_combination)(combination)
                for combination in combination_list
            ]
        )

    (
        return_candidate_nominal_ratio_list,
        return_candidate_composition_list,
        return_input_nominal_ratio_list,
        return_input_compositions_list,
        return_codid_list,
        return_space_group_list,
    ) = zip(*results)

    return_candidate_nominal_ratio_list = list(
        itertools.chain.from_iterable(return_candidate_nominal_ratio_list)
    )

    return_candidate_composition_list = list(
        itertools.chain.from_iterable(return_candidate_composition_list)
    )

    return_input_nominal_ratio_list = list(
        itertools.chain.from_iterable(return_input_nominal_ratio_list)
    )

    return_input_compositions_list = list(
        itertools.chain.from_iterable(return_input_compositions_list)
    )

    return_codid_list = list(itertools.chain.from_iterable(return_codid_list))

    return_space_group_list = list(
        itertools.chain.from_iterable(return_space_group_list)
    )

    return (
        return_candidate_nominal_ratio_list,
        return_candidate_composition_list,
        return_input_nominal_ratio_list,
        return_input_compositions_list,
        return_codid_list,
        return_space_group_list,
    )


def return_chemical_system(cif):
    def get_integer_formula_and_factor(
        self, max_denominator: int = 10000, iupac_ordering: bool = False
    ) -> Tuple[str, float]:
        """
        The default Composition groups together different ox states which is not ideal...
        """
        el_amt = self.as_dict()
        g = gcd_float(list(el_amt.values()), 1 / max_denominator)

        d = {k: round(v / g) for k, v in el_amt.items()}
        (formula, factor) = reduce_formula(d, iupac_ordering=iupac_ordering)
        if formula in Composition.special_formulas:
            formula = Composition.special_formulas[formula]
            factor /= 2
        return formula, factor * g

    # Patched extra functionalities and bug fixes on top of Pymatgen's classes.
    def to_int_dict(self):
        """
        Returns:
            Dict with element symbol and integer amount
        """
        _, factor = self.get_integer_formula_and_factor()
        int_dict = {e: int(a) for e, a in (self / factor).as_dict().items()}
        # be safe: Composition groups together different ox states which is not ideal...
        if not all(
            np.isclose(x * factor, y)
            for x, y in zip(int_dict.values(), self.as_dict().values())
        ):
            raise ValueError("Composition is not rational!")
        return int_dict

    @property
    def inted_composition(self):
        """
        Return Composition instance with integer formula
        """
        _, factor = self.get_integer_formula_and_factor()
        int_comp = self / factor

        # be safe
        int_dict = {e: int(a) for e, a in int_comp.as_dict().items()}
        if not all(
            np.isclose(x * factor, y)
            for x, y in zip(int_dict.values(), self.as_dict().values())
        ):
            raise ValueError("Composition is not rational!")

        return int_comp

    Composition.to_int_dict = to_int_dict
    Composition.get_integer_formula_and_factor = get_integer_formula_and_factor
    Composition.inted_composition = inted_composition

    class PatchedSymmetrizedStructure(SymmetrizedStructure):
        """
        Fixed site_properties display
        """

        def __str__(self):
            outs = [
                "SymmetrizedStructure",
                f"Full Formula ({self.composition.formula})",
                f"Reduced Formula: {self.composition.reduced_formula}",
            ]

            def to_s(x):
                return f"{x:0.6f}"

            outs.append(
                "abc   : "
                + " ".join([to_s(i).rjust(10) for i in self.lattice.abc])
            )
            outs.append(
                "angles: "
                + " ".join([to_s(i).rjust(10) for i in self.lattice.angles])
            )
            if self._charge:
                if self._charge >= 0:
                    outs.append(f"Overall Charge: +{self.charge}")
                else:
                    outs.append(f"Overall Charge: -{self._charge}")
            outs.append(f"Sites ({len(self)})")
            data = []
            props = self.site_properties  # This should be updated!
            keys = sorted(props.keys())
            for i, sites in enumerate(self.equivalent_sites):
                site = sites[0]
                row = [str(i), site.species_string]
                row.extend([to_s(j) for j in site.frac_coords])
                row.append(self.wyckoff_symbols[i])
                for k in keys:
                    row.append(site.properties[k])  # This line
                data.append(row)
            outs.append(
                tabulate(
                    data,
                    headers=["#", "SP", "a", "b", "c", "Wyckoff"] + keys,
                )
            )
            return "\n".join(outs)

    class LabeledStructure(Structure):
        """
        Structure + CIF's _atom_site_label
        """

        def __str__(self):
            return "LabeledStructure\n" + super().__str__()

        @classmethod
        def from_file(  # pylint: disable=arguments-differ
            cls,
            filename,
            primitive=False,
            sort=False,
            merge_tol=0.0,
            symmetrize=False,
        ):
            fname = os.path.basename(filename)
            if not fnmatch(fname.lower(), "*.cif*") and not fnmatch(
                fname.lower(), "*.mcif*"
            ):
                raise ValueError("LabeledStructure only accepts CIFs.")
            instance = super().from_file(
                filename, primitive=primitive, sort=sort, merge_tol=merge_tol
            )

            instance.read_label(filename, symmetrize=symmetrize)
            return instance

        def read_label(
            self,
            cif_filename,
            symprec=const.DEFAULT_SYMPREC,
            angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
            symmetrize=False,
        ):
            """
            Add _atom_site_label as site_properties.
            This is useful for enforcing a certain concentration over
            a group of sites, which may not necessarily consists
            of a single orbit after enlargement to a supercell.

            Args:
                cif_filename (str): Source CIF file.
                symprec (float): Precision for the symmetry operations.

            Raises:
                RuntimeError: If any sites couldn't be matched
                    to one any sites defined within the CIF.
            """
            logging.info(f"Reading _atom_site_label from {cif_filename}")

            def ufloat(string):
                """Remove uncertainties notion from floats (if any)."""
                try:
                    return float(string)
                except ValueError:
                    return float(string.split("(")[0])

            encoding = getattr(io, "LOCALE_ENCODING", "utf8")
            with open(
                cif_filename, "r", encoding=encoding, errors="surrogateescape"
            ) as f:
                parser = CifParser.from_string(f.read())

            # Since Structure only takes the first structure inside a CIF, do the same.
            cif_dict = list(parser.as_dict().values())[0]
            labels = cif_dict["_atom_site_label"]
            x_list = map(ufloat, cif_dict["_atom_site_fract_x"])
            y_list = map(ufloat, cif_dict["_atom_site_fract_y"])
            z_list = map(ufloat, cif_dict["_atom_site_fract_z"])
            coords = [(x, y, z) for x, y, z in zip(x_list, y_list, z_list)]

            # Merge labels to allow multiple references.
            cif_sites = []
            for coord, zipgroup in itertools.groupby(
                zip(coords, labels), key=lambda x: x[0]
            ):
                labels = tuple(sorted({x[1] for x in zipgroup}))
                cif_sites.append(
                    PeriodicSite(
                        "X",
                        coord,
                        self.lattice,
                        properties={"_atom_site_label": labels},
                    )
                )

            # Find equivalent sites.
            if symmetrize:
                symm_ops = SpacegroupAnalyzer(
                    self, symprec, angle_tolerance
                ).get_space_group_operations()
            else:
                # A bit of trick.
                parser.data = cif_dict
                # Spacegroup symbol and number are not important here.
                symm_ops = SpacegroupOperations(
                    0, 0, parser.get_symops(parser)
                )

            coords = [x.frac_coords for x in self.sites]
            cif_coords = [x.frac_coords for x in cif_sites]

            # List of coordinates that are equivalent to this site
            o_cif_coords = [
                symmop.operate_multi(cif_coords) for symmop in symm_ops
            ]
            o_cif_coords = np.swapaxes(np.stack(o_cif_coords), 0, 1)
            o_cif_coords = [
                np.unique(np.mod(x, 1), axis=0) for x in o_cif_coords
            ]

            for site in tqdm(
                self.sites,
                desc="Matching CIF labels",
                **const.TQDM_CONF,
                disable=True,
            ):
                equivalent = [
                    in_coord_list_pbc(o, site.frac_coords)
                    for o in o_cif_coords
                ]

                try:
                    equivalent_site = cif_sites[equivalent.index(True)]
                    site.properties[
                        "_atom_site_label"
                    ] = equivalent_site.properties["_atom_site_label"]
                except ValueError as exc:
                    raise RuntimeError("CIF-Structure mismatch.") from exc

    class PatchedSpacegroupAnalyzer(SpacegroupAnalyzer):
        """
        Patched to use PatchedSymmetrizedStructure in get_symmetrized_structure()
        """

        def get_symmetrized_structure(self):
            """
            Use PatchedSymmetrizedStructure
            """
            ds = self.get_symmetry_dataset()
            sg = SpacegroupOperations(
                self.get_space_group_symbol(),
                self.get_space_group_number(),
                self.get_symmetry_operations(),
            )
            return PatchedSymmetrizedStructure(
                self._structure, sg, ds["equivalent_atoms"], ds["wyckoffs"]
            )

    def beautiful_inted_formula(composition):
        inted_element_composition = (
            composition.inted_composition.element_composition
        )
        # Because buggy formula!
        sym_amt = inted_element_composition.get_el_amt_dict()
        syms = sorted(sym_amt.keys(), key=lambda sym: get_el_sp(sym).X)
        formula = [
            s + formula_double_format(int(sym_amt[s]), False) for s in syms
        ]
        return "".join(formula)

    cifid, _ = os.path.splitext(os.path.basename(cif))

    try:
        parser = CifParser(cif)

        key = list(parser.as_dict().keys())[0]
        cif_dict = parser.as_dict()[key]
        chem_formula = cif_dict["_chemical_formula_sum"].replace(" ", "")

        structure = LabeledStructure.from_file(cif)
        sym_structure = PatchedSpacegroupAnalyzer(
            structure,
            symprec=const.DEFAULT_SYMPREC,
            angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        ).get_symmetrized_structure()
        sga = PatchedSpacegroupAnalyzer(
            structure,
            symprec=SYM_TOLERANCE,
            angle_tolerance=ANGLE_TOLERANCE,
        )
        if "_space_group_IT_number" in cif_dict:
            cif_sg_number = cif_dict["_space_group_IT_number"]
        else:
            cif_sg_number = cif_dict["_symmetry_Int_Tables_number"]
        if sga.get_space_group_number() != int(cif_sg_number):
            raise ValueError(
                "the SG in the cif file is not consistent with the one detected by spglib."
            )

        lattice_type = sga.get_lattice_type()
        point_group = sga.get_point_group_symbol()
        space_group_num = sga.get_space_group_number()
        space_group = sga.get_space_group_symbol()

        sub_chem_formula = beautiful_inted_formula(sym_structure.composition)

        chemical_system = sym_structure.composition.chemical_system
        element_key = hashlib.md5(chemical_system.encode("utf-8")).hexdigest()
        element_key = chemical_system

        record = [
            cifid,
            lattice_type,
            point_group,
            space_group_num,
            space_group,
            chem_formula,
            sub_chem_formula,
        ]

        return element_key, record

    except (ValueError, TypeError, KeyError) as e:
        e
        # print(cifid)
        # print(e)
        return ("NA", "NA")

    """
    except Exception as e:
        print(f"cifid={cifid}, {e}")
        return ("NA", "NA")
    """
