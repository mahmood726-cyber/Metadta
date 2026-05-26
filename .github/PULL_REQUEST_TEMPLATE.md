## What this PR changes

<!-- One paragraph: the change, not the motivation. -->

## Why

<!-- The bug, gap, or numerical discrepancy this fixes. Include DOI /
issue link if applicable. -->

## How I verified

- [ ] `python -m pytest -q` passes locally
- [ ] `ruff check .` passes (pre-commit will run this on commit)
- [ ] Numerical results (if any) verified against R reference or
      previous tag — paste the diff here:

```
# pasted comparison
```

## Checklist

- [ ] No hardcoded local paths in committed code
      (`C:\Users\...`, `/home/<user>/...`)
- [ ] If new dependency added, pinned in `requirements.txt` or
      `pyproject.toml`
- [ ] If user-facing string changed, dashboard render still looks right
- [ ] Commit messages explain the *why*, not just the *what*
