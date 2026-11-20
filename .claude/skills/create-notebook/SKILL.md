---
name: create-notebook
description: End-to-end check of Providentia's Library mode. Asks the user (pop-up questions) for the conf, section, subsection, plots and plot options, builds a notebook with the standard setup (import providentia, Providentia(conf, section, subsection), print_config, load) plus one cell per plot in notebooks/, runs it if the user says so, and gives its path. Use when the user wants to try a library feature, plot or fix in a notebook.
---

# Create a Library mode notebook

Every notebook is: import providentia, create the class, print_config(), load(), then
the test cells. The setup is written by notebook.py, never by hand. Ignore anything typed
after /create-notebook: always ask the questions below, with the AskUserQuestion tool.

## 0. Environment
Take the environment name from the local-run line of the launcher, never from memory, as
it changes between branches:
`ENV_NAME=$(grep -m1 '^ENV_NAME=' bin/providentia | cut -d'"' -f2)`, then
`ENV=$(dirname $(conda run -n $ENV_NAME which python))`. If either is empty or conda cannot
find the environment, stop and tell the user. Never use the shell's `python` (base
anaconda can have an older matplotlib).

## 1. Ask, in rounds (each round needs the answers of the one before)
Tab titles (`header`, 12 characters at most) are short title-case names: "Config",
"Section", "Subsection", "Plots", and for plot options the plot name ("Heatmap",
"Statsummary"). Option descriptions are left empty when the label says it all.

### Round 1: Configuration
"Which configuration? (last 3 edited)" The 3 most recently modified
`configurations/*.conf` (`ls -t`), shown without `.conf`; the user types any other name
in Other. Stop and ask again if it does not exist.

### Round 2: Section, subsection and plots
Read from that conf (`[section]`, `[[subsection]]`, the section's `plots =` line):
- "Which section?" "None" first (the library then takes the first section), then up to
  3 sections; others in Other. If the conf has only one section, do not ask: use that
  section and pass it with `--section`.
- "Which subsection?" "None" first (the library then takes the first subsection), then
  up to 3 subsections of the chosen section, or of the first section if the answer was
  "None".
- "Which plots?" multiSelect: the plot types in the section's `plots =` line first,
  then common ones (timeseries, distribution, map, scatter) up to 4. Plot types are the
  ones in `docs/source/Library.md`.

### Round 3: Plot options
One multiSelect question per chosen plot (4 per pop-up, more rounds if needed):
- "Options for <plot>?" "None" first, then the first 3 options of that plot type in
  `settings/plot_characteristics.yaml`, under `<plot>` → `library` → `plot_options`
  (`<plot>` without the statistic: `heatmap` for `heatmap-Mean`; if there is no `library`
  block, the plot type's own `plot_options`). The rest are typed in Other. Skip the
  question for plot types with no options.

### Round 4: Run
Always in its own pop-up after the plot options: "Run it now?" (tab title "Run"),
options "Yes" (the plots are in the notebook when it is opened) and "No". The pop-up
always shows Other; if the answer typed there is not clearly yes or no, ask again.

### Answers typed in Other
In multiSelect questions, Other can hold several items separated by commas; the question
text ends with "(Other: several, comma-separated)". Split it and add each item to the
ticked options, checking plot types against `docs/source/Library.md` and plot options
against the list in `settings/plot_characteristics.yaml` used in Round 3. Keep the valid
items; for any unknown one, ask the same pop-up question again, with the question text
naming the unknown item and listing every valid plot type (or every option of that plot
type).
Anything else the user types in Other (statistics, filters, set_config, free code) goes
into the test cells too; check it against `providentia/library.py` and
`docs/source/Library.md`, and ask if it is not clear.

## 2. Write the test cells
In a temporary `.py` file, cells separated by `# %%`, one cell per plot:
`provi.plot('<plot>')`, or `provi.plot('<plot>', plot_options=['annotate', 'bias'])`.

## 3. Build
`$ENV/python .claude/skills/create-notebook/notebook.py build --conf <conf>.conf [--section <s>]
[--subsection <ss>] --cells <cells.py> --out notebooks/<YYYYMMDD>_<conf>.ipynb`
where `<conf>` is the configuration name without `.conf` (e.g. `20261111_debug.ipynb`).
Leave out `--section` / `--subsection` when the answer was "None". It prints the notebook
path, with `_2`, `_3` added if the name exists; an existing notebook is never replaced.

## 4. Run, if the answer was "Yes"
From the repo root, with `PYTHONPATH` set to it as `bin/providentia --notebook` does (the
kernel starts in `notebooks/`, where `import providentia` fails without it):
`PYTHONPATH=$(pwd) $ENV/jupyter nbconvert --to notebook --execute --inplace --allow-errors
--ExecutePreprocessor.timeout=1800 <notebook>`
If it complains that folders in the home directory are not writable (Claude's sandbox),
set MPLCONFIGDIR, XDG_CACHE_HOME and IPYTHONDIR to a temporary folder. Then
`grep -c '"output_type": "error"' <notebook>` tells whether a cell raised an error.

## 5. Report
Tell the user the notebook path and its cells, and whether it was run (and if a cell
raised an error, which one). Open it with `./bin/providentia --notebook` (or in a
JupyterLab that is already open). It is in `notebooks/`, which git tracks, so it shows in
`git status` until it is deleted or committed. A run rewrites the filetree JSONs in
`settings/internal/`; that is normal and not to be committed.

Notes:
- Do not edit files in `providentia/` from this skill. If the user reports a bug from the
  notebook, report it back and ask.
