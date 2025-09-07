# Contributing

1. Create a virtual environment and install in editable mode:

   python -m venv .venv
   source .venv/bin/activate
   pip install -e .[dev]

2. Run the CLI:

   vamos-precision map
   vamos-precision cluster
   vamos-precision gsea-prep
   vamos-precision gsea-run
   vamos-precision ml
   vamos-precision logodds
   vamos-precision prioritize
   vamos-precision drug-id
   vamos-precision permute

3. Tests:

   pytest -q
