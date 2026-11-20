# Skills

Skills are step-by-step instructions, kept in the repository, that an AI coding assistant can follow to do a common Providentia task for you. Each skill is a folder in `.claude/skills/` with a `SKILL.md` file. To use one, type its name after a slash in your assistant, for example `/create-notebook`. Skills have only been tried with one assistant so far: others may look for them in a different folder, or ask their questions differently.

## Create a Library mode notebook

The `/create-notebook` skill saves writing a test notebook by hand. Type `/create-notebook` and answer the questions that pop up: the configuration file, the section (not asked when the file has only one) and subsection (or "None", to take the first one), the plots and their plot options.

It creates a notebook in `notebooks/` that imports Providentia, creates the class with your configuration, loads the data and then makes the plots you chose. The last question asks whether to run it now: if you say yes, it is run in a new Python process, so the plots are already there when you open it with `./bin/providentia --notebook` (see [Notebooks](Notebooks.md)). Delete it when you are done, as git tracks the `notebooks/` folder.
