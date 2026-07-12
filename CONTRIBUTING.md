# Contributing

Contributions are welcome through GitHub issues and pull requests.

## Development setup

```bash
git clone https://github.com/andrestrocyte/pymutscan.git
cd pymutscan
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,notebooks]"
python scripts/generate_synthetic_mapseq.py
python -m unittest discover -s tests -v
ruff check src tests scripts
```

Run the notebooks exactly as CI does:

```bash
for notebook in notebooks/*.ipynb; do
  jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=180 "$notebook"
done
```

## Pull requests

- Add tests for behavioral changes.
- Preserve raw-sequence provenance and mapping auditability.
- Document any change to clustering semantics or tie-breaking.
- Do not commit identifiable or unpublished experimental sequencing data.
- Keep example data deterministic and synthetic.

