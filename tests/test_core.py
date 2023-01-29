# Copyright (c) cifcomb Development Team.
# Distributed under the terms of the MIT License.

# python modules
from helper import chdir
from pymatgen.core import Composition

from cifcomb.core import decompose_composition


def test_decompose_composition():
    target_composition = Composition("CaTiO3")
    input_composition_list = [Composition(i) for i in ["CaO", "TiO2"]]
    max_multiplication = 10
    target_nominal_ratio, input_nominal_ratio_list = decompose_composition(
        target_composition=target_composition,
        input_composition_list=input_composition_list,
        max_multiplication=max_multiplication,
    )

    assert target_nominal_ratio == 1
    assert input_nominal_ratio_list == [1, 1]


@chdir("../examples")
def test_example2():
    # to be implemented
    assert True


if __name__ == "__main__":
    test_decompose_composition()
