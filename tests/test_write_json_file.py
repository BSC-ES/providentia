import json
import os

import pytest

from providentia.read_aux import write_json_file


def test_write_json_file_creates_file(tmp_path):
    path = str(tmp_path / "filetree.json")

    write_json_file(path, {"EBAS": ["sconco3"]})

    with open(path) as json_file:
        assert json.load(json_file) == {"EBAS": ["sconco3"]}


def test_write_json_file_replaces_file_without_leaving_temporary_file(tmp_path):
    path = str(tmp_path / "filetree.json")

    write_json_file(path, {"old": True})
    write_json_file(path, {"new": True})

    with open(path) as json_file:
        assert json.load(json_file) == {"new": True}
    assert os.listdir(tmp_path) == ["filetree.json"]


def test_write_json_file_keeps_old_file_if_writing_fails(tmp_path):
    path = str(tmp_path / "filetree.json")
    write_json_file(path, {"old": True})

    # sets cannot be written to JSON, so writing stops halfway
    with pytest.raises(TypeError):
        write_json_file(path, {"new": True, "invalid": {1, 2}})

    with open(path) as json_file:
        assert json.load(json_file) == {"old": True}
