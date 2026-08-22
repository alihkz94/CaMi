#!/usr/bin/env python3
"""
=============================================================================
check_pins.py — the two statistics environments must describe the same software
-----------------------------------------------------------------------------
CaMi installs the scientific Python stack twice, by two different routes:

    env/requirements-stats.txt   pip, used by the container and by the venv
                                 that scripts/analysis/setup_stats_env.sh builds
    env/stats-conda.yml          conda + a pip section, used by -profile conda

Both exist for reasons written down in those files. The danger is not that they
exist; it is that someone bumps scipy in one of them and not the other, and then
two runs that both claim to be CaMi 1.0.0 produce two different p-values with
nothing to show why.

This script fails the build when the pinned versions disagree. It is deliberately
strict: it accepts only `pkg==version`, because a range in either file would make
the question unanswerable rather than merely wrong.

USE:
    python3 scripts/containers/check_pins.py
Exit status 0 when they agree, 1 when they do not.
=============================================================================
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
REQ = REPO / "env" / "requirements-stats.txt"
YML = REPO / "env" / "stats-conda.yml"

PIN = re.compile(r"^([A-Za-z0-9_.\-]+)==([A-Za-z0-9_.\-+!]+)$")


def parse_requirements(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        m = PIN.match(line)
        if not m:
            sys.exit(f"{path.name}:{lineno}: not an exact pin: {line!r}")
        out[m.group(1).lower()] = m.group(2)
    return out


def parse_conda_pip_section(path: Path) -> dict[str, str]:
    """
    Read only the `pip:` list. A full YAML parse would need PyYAML, and this file
    is checked by CI on a bare runner where adding a dependency to read a
    dependency file is the wrong trade.
    """
    out: dict[str, str] = {}
    in_pip = False
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if re.match(r"^\s*-\s*pip:\s*$", line):
            in_pip = True
            continue
        if in_pip:
            m = re.match(r"^\s+-\s+(.+?)\s*$", line)
            if not m:
                in_pip = False          # the pip list ended
                continue
            pin = PIN.match(m.group(1))
            if not pin:
                sys.exit(f"{path.name}:{lineno}: not an exact pin: {m.group(1)!r}")
            out[pin.group(1).lower()] = pin.group(2)
    return out


def main() -> int:
    req = parse_requirements(REQ)
    yml = parse_conda_pip_section(YML)

    if not req:
        sys.exit(f"{REQ} lists no packages at all — that cannot be right.")
    if not yml:
        sys.exit(f"{YML} has no pip: section — that cannot be right.")

    problems: list[str] = []
    for pkg in sorted(set(req) | set(yml)):
        a, b = req.get(pkg), yml.get(pkg)
        if a is None:
            problems.append(f"  {pkg}: only in {YML.name} ({b})")
        elif b is None:
            problems.append(f"  {pkg}: only in {REQ.name} ({a})")
        elif a != b:
            problems.append(f"  {pkg}: {REQ.name} says {a}, {YML.name} says {b}")

    if problems:
        print("The two statistics environments disagree:\n")
        print("\n".join(problems))
        print(
            "\nFix BOTH files. A container and a conda profile that claim to be the"
            "\nsame version of CaMi must not install different software."
        )
        return 1

    print(f"OK — {len(req)} pinned package(s) agree between {REQ.name} and {YML.name}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
