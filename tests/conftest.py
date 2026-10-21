import os

# make get_machine() return "github" whenever the test suite runs,
# so tests use the same paths locally as in GitHub Actions
os.environ.setdefault("GITHUB_ACTIONS", "true")
