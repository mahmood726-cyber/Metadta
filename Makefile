# Convenience targets. `make help` lists them.
.PHONY: help install test cov lint fmt clean pre-commit-install pre-commit-run

PY ?= python

help:	## Show this help message
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install:	## Install runtime + test deps + pre-commit
	$(PY) -m pip install --upgrade pip
	if [ -f requirements.txt ]; then $(PY) -m pip install -r requirements.txt; fi
	if [ -f pyproject.toml ]; then $(PY) -m pip install -e ".[dev]" 2>/dev/null || $(PY) -m pip install -e .; fi
	$(PY) -m pip install ruff pytest pytest-cov pre-commit
	pre-commit install

test:	## Run pytest (terse)
	$(PY) -m pytest -q --tb=short

cov:	## Run pytest with coverage, HTML report into htmlcov/
	$(PY) -m pytest -q --cov=. --cov-report=term-missing --cov-report=html
	@echo "Open htmlcov/index.html for the coverage report"

lint:	## Run ruff check (no fixes)
	$(PY) -m ruff check .

fmt:	## Run ruff format + auto-fixable lint
	$(PY) -m ruff format .
	$(PY) -m ruff check --fix .

clean:	## Remove caches and coverage artifacts
	rm -rf .pytest_cache .ruff_cache .mypy_cache htmlcov coverage.xml .coverage
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name "*.egg-info" -prune -exec rm -rf {} +

pre-commit-install:	## Install pre-commit hooks
	pre-commit install

pre-commit-run:	## Run all pre-commit hooks against all files
	pre-commit run --all-files
