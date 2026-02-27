# HUMAnN4 Tools - Agent Guidelines

## Build/Test/Lint Commands
- **Build:** `pip install -e .` or `pip install -e ".[dev]"` (with dev dependencies)
- **Environment:** `python conda_setup.py [--name ENV_NAME] [--python PYTHON_VERSION]`
- **Lint:** `black .` (line length 100), `isort .`, `flake8 .`, `mypy .`
- **Test:** `pytest` or `pytest tests/`, `pytest --cov=humann4_tools`

## Code Style
- **Line Length:** 100 characters
- **Formatting:** Black with Python 3.8+ targets
- **Type Hints:** Used but not strictly enforced; mypy configured with `warn_return_any = true`
- **Imports:** Sorted with isort using Black compatibility profile
- **Error Handling:** Use comprehensive logging with `logger.py` functions
- **Naming Convention:** snake_case for functions/variables, CamelCase for classes
- **Documentation:** Clear docstrings with Parameters and Returns sections
- **Command Structure:** Follow existing CLI patterns with argparse

## Project Structure
- Core package: `humann4_tools/`
- Main workflows in: `preprocessing/`, `analysis/`, `humann4/` 
- Utils in: `utils/` (file operations, metadata handling, resources)
- Entry points defined in pyproject.toml under `[project.scripts]`