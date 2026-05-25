# Security policy

## Supported versions

| Version | Supported |
| --- | --- |
| Latest commit on `main` / `master` | yes |
| Tagged releases (if any) | most recent only |
| Older commits | no |

This is a research-software repository — it does not ship a long-term
support track. Treat the tip of the default branch as canonical.

## Reporting a vulnerability

If you find a vulnerability in this code, please **do not** open a
public issue. Instead:

1. Email **mahmood.ahmad2@nhs.net** with subject prefix `[SECURITY]`,
   or
2. Use GitHub's [private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories/guidance-on-reporting-and-writing/privately-reporting-a-security-vulnerability)
   feature on this repository (Security tab → "Report a vulnerability").

I aim to acknowledge security reports within 7 days and to ship a fix
or mitigation within 30 days for confirmed issues. Co-disclosure timing
is negotiable.

## Scope

This repository is a browser-based / Python research-software tool.
In-scope vulnerabilities include:

- Code-injection / XSS in the rendered dashboard
- Unsafe deserialization of user-supplied data (JSON / CSV / YAML)
- Hardcoded secrets or credentials
- Path-traversal in any file-read path
- Dependency CVEs that affect this repo's runtime

Out of scope: vulnerabilities in third-party services this tool talks
to (CT.gov, CrossRef, PubMed) — report those upstream.

## Automated checks

This repository runs:

- **CodeQL** (`.github/workflows/codeql.yml`) — weekly + on every push/PR
- **pip-audit** (manual; declared deps verified at commit time)
- **Dependabot** (`.github/dependabot.yml`) — weekly dep update PRs
- **Sentinel pre-push** (per the upstream portfolio's pre-push hooks)

See those workflows for the current rule sets.
