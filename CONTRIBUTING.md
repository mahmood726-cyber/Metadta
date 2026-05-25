# Contributing

Thanks for your interest in this project. This is part of the
[E156 micro-paper portfolio](https://github.com/mahmood726-cyber/e156) —
a collection of browser-based, dependency-light research-software tools.

## Quick start

```bash
git clone <repo-url>
cd <repo-name>

# Install dev deps
pip install -r requirements.txt 2>/dev/null || pip install pytest pytest-cov

# Install pre-commit hooks (ruff lint + format on every commit)
pip install pre-commit
pre-commit install

# Run the test suite
python -m pytest -q
```

## What kinds of contributions are welcome

- **Bug reports** with a minimal reproducer
- **Numerical-parity findings** — if results disagree with R metafor / meta /
  mada or similar reference implementations
- **Documentation fixes** — typos, broken links, unclear instructions
- **Test coverage** for currently-untested branches
- **Translations** of dashboards into additional languages (the existing
  dashboards aim to be language-portable)

## What's out of scope (most of the time)

- Reformatting / restyling pull requests that don't fix a bug or add coverage
  (the pre-commit + ruff config already handles formatting consistency)
- Adding heavy runtime dependencies (the dashboards aim to stay
  dependency-light and browser-runnable)

## Pull-request checklist

- [ ] `python -m pytest -q` passes locally
- [ ] `ruff check .` passes (pre-commit will run this automatically)
- [ ] If you added a Python dependency, it's pinned in `requirements.txt`
      (or `pyproject.toml`)
- [ ] Commit messages explain the *why*, not just the *what*
- [ ] No hardcoded local paths (`C:\Users\...`, `/home/<user>/...`) in
      committed code

## Authorship for E156 micro-paper submissions

If this repo backs a published or submission-ready E156 micro-paper, the
authorship convention is documented in `E156-PROTOCOL.md`. In short:
the student rewriter is first author; their faculty supervisor is
last/senior author; Mahmood Ahmad is middle author only.

## Reporting security issues

See [`SECURITY.md`](SECURITY.md).

## License

By contributing, you agree your contributions are licensed under the
same MIT License that covers this repository (see `LICENSE`).
