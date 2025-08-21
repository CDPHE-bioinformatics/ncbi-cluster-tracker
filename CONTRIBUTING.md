We welcome contributions through issues and pull requests!

# Poetry instructions

[Install](https://python-poetry.org/docs/#installation) the `poetry` command-line tool for packaging and dependency management.

Ensure you have Python >= 3.12 installed in your environment. If `poetry` is unable to find the correct Python version, you can manually set with `poetry env use <path to python3.12>`.

From the base directory of the cloned/forked repository, install dependencies for `ncbi-cluster-tracker` with
```
poetry install
```

Run the tool with local changes using
```
poetry run ncbi-cluster-tracker ...
```

To create and publish a distributable package to PyPI, ensure the `version` under `[project]` is updated in the `pyproject.toml` and run

```
poetry build
poetry publish
```

# Testing
Example dataset:

```
poetry run ncbi-cluster-tracker --keep-snp-files --retry --out-dir tests/data/outputs/20250416_180514 tests/data/sample_sheet.csv
```

Unit tests:

```
python -m unittest tests/tests.py
```

Static type checking with [mypy](https://mypy.readthedocs.io/en/stable/):

```
mypy ncbi_cluster_tracker --strict
```