# AGENTS.md

## Overview

This is a simple pipeline to download the constituent datasets of the Unified
Genealogy (https://github.com/awohns/unified_genealogy_paper,
https://www.science.org/doi/10.1126/science.abi8264) focusing on just SGDP and
1000 Genomes Project for now, and run tsinfer 1.0 with both sources. It has to
be done in GRCh37 because no phased GRCh38 SGDP dataset is currently available.
We also need to download the ancestral allele data and output it as a zarr and
recombination maps too. The Snakemake pipeline just has to download the data for
chromosomes specified in its config.yaml from the source URLs specified there,
then run bio2zarr to convert them into the right format for tsinfer 1.0. For now
(and possibly forever) don't integrate tsinfer into the Snakemake pipeline, it
just has to assemble everything needed.

The snakemake architecture should take inspiration from
https://github.com/benjeffery/tsinfer-snakemake/, although bear in mind this is
for an old version of tsinfer. The organisation of a single Snakefile and
steps.py with functions it calls should be in place here too. 

## Git usage

- Don't include Co-authored-By lines in git commits.

## Code Style

- Do not be overly defensive - defend only against circumstances that can
  occur within the current codebase.
- Do not make production code more complex for the sake of minimising 
  changes to the test suite. Simplicity and clarity of the production code 
  is imperative.
- Do not combine multiple complex operations in a single statement. Prefer
  to keep a single operation per statement, and use intermediate variables
  as a form of documentation. For example:
  ```python
  # Bad — multiple operations in one expression
  result = sorted(k for k, v in mapping.items() if v in set(x.name for x in sources))

  # Good — intermediate variable makes intent clear
  source_names = {x.name for x in sources}
  result = sorted(k for k, v in mapping.items() if v in source_names)
  ```
- Prefer dataclasses over tuples when returning multiple values.
- Use explicit `None` comparisons: `if x is not None` not `if x`.
- Import all modules at the top of the file, not inside functions or methods.
- Prefer importing a module and using module.function instead of
  using ``from module import function``. This applies to intra-package
  imports too: use ``from . import config`` then ``config.X``, not
  ``from .config import X``. Exceptions: ``from typing import ...`` is
  acceptable; ``from .X import Y`` is acceptable in ``__init__.py`` for
  defining the public API. Use ``import concurrent.futures as cf``.
- Use idiomatic pathlib.Path operations instead of os.path operations.
- When a parameter has a computed default derived from another parameter,
  compute it once at the point of use (the leaf function), not at every
  layer in the call chain. Pass `None` through intermediate layers.
- Zarr v3 is used (dependency: `zarr>=3`). Do not use Zarr v2 APIs.
- Use PEP 604 union syntax: `int | None`, not `Optional[int]`.
- One `logger = logging.getLogger(__name__)` per module at top level.

## Tool use

- Use `uv run` for all Python tooling (never bare `python -m`)

## Testing

- Organise tests in classes, not flat functions. Use pytest fixtures for setup.
- Test helpers are in `tests/helpers.py` (e.g., `make_sample_vcz`, `make_ancestor_vcz`)
- `tests/algorithm.py` contains Python reference implementations used to verify C code
- `msprime` is used to simulate test data
- Run the test suite with coverage before committing to ensure that new code is
  fully covered by tests.


## Coding style and quality

- Use discernment and good judgement to write code that is correct, maintainable, and efficient, while avoiding unnecessary complexity or over-engineering.
- Use underscores with lower case for variable names in CSVs and for files and folders except .agents, where you can use the - separation convention.
- Prioritise semantic correctness, well-defined invariants, and clear error handling over micro-optimizations or clever tricks.
- Maintain inter-linked docstrings and use comments where they add context and explanation, but avoid redundant or obvious comments that do not enhance understanding.
- Structure code to be modular, with coherent responsibilities, good separation of concerns, and an architecture that can be extended without large-scale rewrites.
- Apply DRY carefully: factor out shared logic when it improves clarity and maintenance, but avoid over-abstraction or generic frameworks that obscure intent.
- Aim for readable, explicit code with simple control flow and minimal “magic”; prefer straightforward implementations to compact but opaque idioms.
- Consider algorithmic time and memory complexity when choosing data structures and approaches, especially for code on critical paths or large problem sizes.
- Keep implementations concise but not cryptic; remove dead or unused code and avoid unnecessary indirection.
- For tests, prefer simple, direct assertions that fail fast; test code should not be defensive, over-general, or attempt to mask or recover from failures.
- For APIs and function interfaces prefer defined contracts with clear error states over flexible but loosely defined behavior that attempts to infer intent or handle edge cases implicitly.
- Where there is likely to be a long running process and it can be parallelised use multiprocessing with queues to distribute work across all available CPU cores.
- Don't bend over backwards to make anything backwards compatible with previous implementations. 
