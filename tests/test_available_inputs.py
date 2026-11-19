import importlib.util
from os.path import join

import yaml

from providentia.auxiliar import CURRENT_PATH, PROVIDENTIA_ROOT


def test_available_species_in_ghost_standards():
    """Every species of available_species_per_matrix exists in the GHOST 1.6 standards, under the matrix it is listed in."""

    # load the GHOST 1.6 standards by path, so another version imported earlier is not used
    spec = importlib.util.spec_from_file_location(
        "GHOST_standards_1_6",
        join(CURRENT_PATH, "dependencies", "GHOST_standards", "1.6", "GHOST_standards.py"),
    )
    ghost_standards = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ghost_standards)
    ghost_matrices = {
        param_dict["bsc_parameter_name"]: param_dict["matrix"]
        for param_dict in ghost_standards.standard_parameters.values()
    }

    with open(join(PROVIDENTIA_ROOT, "settings", "available_inputs.yaml")) as f:
        available_species_per_matrix = yaml.safe_load(f)["available_species_per_matrix"]

    wrong = [
        f"{speci} (listed under {matrix}, GHOST 1.6: {ghost_matrices.get(speci, 'not found')})"
        for matrix, species in available_species_per_matrix.items()
        for speci in species
        if ghost_matrices.get(speci) != matrix
    ]
    assert not wrong, f"Species not in GHOST 1.6 under their matrix: {', '.join(wrong)}"
