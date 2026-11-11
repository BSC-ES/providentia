"""Build the notebooks of the /create-notebook skill: the standard Library mode setup cells
(import, class, print_config, load), then the test cells.

Run it with the python of the environment named by ENV_NAME in bin/providentia (local
runs), which has nbformat.
"""

import argparse
import os
import re
import sys

import nbformat

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))

SETUP_CELLS = ["# import Providentia\nimport providentia as prv",
               "# magic to make plots appear inline\n%matplotlib inline"]


def unique_path(path):
    """Return path, or path with _2, _3... before the extension if it already exists."""
    base, ext = os.path.splitext(path)
    n = 2
    while os.path.exists(path):
        path = f"{base}_{n}{ext}"
        n += 1
    return path


def build(args):
    """Write the notebook: setup cells, then the cells of args.cells split on '# %%'."""
    conf = args.conf
    if not os.path.exists(conf):
        conf = os.path.join(REPO_ROOT, "configurations", args.conf)
    if not os.path.exists(conf):
        sys.exit(f"Configuration not found: {args.conf} (also looked in configurations/)")
    conf = os.path.abspath(conf)
    # like the example notebooks, a conf in configurations/ is passed by name only
    if os.path.dirname(conf) == os.path.join(REPO_ROOT, "configurations"):
        conf = os.path.basename(conf)
    init_args = [repr(conf)]
    if args.section:
        init_args.append(f"section={args.section!r}")
    if args.subsection:
        init_args.append(f"subsection={args.subsection!r}")
    init = (f"# initialise objects\nprovi = prv.Providentia({', '.join(init_args)})\n"
            "provi.print_config()")

    with open(args.cells) as f:
        test_cells = [c.strip() for c in re.split(r"^# %%.*$", f.read(), flags=re.M)]

    nb = nbformat.v4.new_notebook()
    nb.cells = [nbformat.v4.new_code_cell(c) for c in SETUP_CELLS]
    nb.cells += [nbformat.v4.new_code_cell(init),
                 nbformat.v4.new_code_cell("provi.load()")]
    nb.cells += [nbformat.v4.new_code_cell(c) for c in test_cells if c]

    out = unique_path(os.path.abspath(args.out))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "x") as f:
        nbformat.write(nb, f)
    print(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_build = sub.add_parser("build", help="write a notebook with the standard setup")
    p_build.add_argument("--conf", required=True, help="path to the .conf file")
    p_build.add_argument("--section", help="omit to let the library take the first section")
    p_build.add_argument("--subsection")
    p_build.add_argument("--cells", required=True, help="file with the test cells, split on '# %%'")
    p_build.add_argument("--out", required=True, help="notebook path; _2, _3... added if it exists")
    p_build.set_defaults(func=build)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
