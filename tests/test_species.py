import copy
import math
from os.path import join
from types import SimpleNamespace

import yaml

from providentia.auxiliar import PROVIDENTIA_ROOT, load_ghost_species, load_species_reference
from providentia.read_aux import apply_ghost_version_fields

FIELDS = [
    "long_parameter_name",
    "matrix",
    "standard_units",
    "units_quantity",
    "chemical_formula",
    "chemical_formula_charge",
    "extreme_lower_limit",
    "extreme_upper_limit",
    "minimum_resolution",
]


def load_reference():
    """Load the species reference of Providentia (settings/species.yaml)"""

    with open(join(PROVIDENTIA_ROOT, "settings", "species.yaml")) as f:
        return yaml.safe_load(f)


def same_value(first, second):
    """Whether two values are equal, taking two NaN as equal"""

    if isinstance(first, float) and isinstance(second, float):
        if math.isnan(first) and math.isnan(second):
            return True
    return first == second


def test_every_ghost_species_is_in_the_reference():
    missing = sorted(set(load_ghost_species("1.6")) - set(load_reference()))
    assert not missing, f"GHOST 1.6 species missing in settings/species.yaml: {', '.join(missing)}"


def test_ghost_species_of_the_reference_match_ghost():
    """If GHOST 1.6 changes, the GHOST species of settings/species.yaml have to be regenerated"""

    ghost_species = load_ghost_species("1.6")
    reference = load_reference()

    different = [
        f"{speci}.{field} (reference: {reference[speci][field]!r}, GHOST 1.6: {ghost_species[speci][field]!r})"
        for speci in sorted(set(reference) & set(ghost_species))
        for field in FIELDS
        if not same_value(reference[speci][field], ghost_species[speci][field])
    ]
    assert not different, f"settings/species.yaml differs from GHOST 1.6: {'; '.join(different)}"


def test_every_species_of_the_reference_has_its_fields():
    """Also species that are not in GHOST, added by hand"""

    wrong = [
        speci
        for speci, fields in load_reference().items()
        if (sorted(fields) != sorted(FIELDS)) or not fields["matrix"]
    ]
    assert not wrong, f"species without exactly the fields {FIELDS}, or without a matrix: {', '.join(wrong)}"


def test_ghost_observations_take_the_fields_of_their_version_but_the_matrix():
    """precmsa is in mg l-1 in GHOST 1.5 and in ug l-1 in the reference; sphaero is in the asy
    matrix in GHOST 1.5 and in size in the reference"""

    species_reference = load_species_reference()
    instance = SimpleNamespace(
        species_reference=species_reference,
        parameter_dictionary=copy.deepcopy(species_reference),
        ghost_parameter_dictionary=load_ghost_species("1.5"),
        reading_ghost=True,
    )

    apply_ghost_version_fields(instance)
    assert instance.parameter_dictionary["precmsa"]["standard_units"] == "mg l-1"
    assert instance.parameter_dictionary["sphaero"]["matrix"] == "size"
    # species that GHOST 1.5 does not have keep the reference's fields
    assert instance.parameter_dictionary["sconcald2"] == species_reference["sconcald2"]

    instance.reading_ghost = False
    apply_ghost_version_fields(instance)
    assert instance.parameter_dictionary["precmsa"]["standard_units"] == "ug l-1"
    assert instance.parameter_dictionary["sphaero"]["matrix"] == "size"
