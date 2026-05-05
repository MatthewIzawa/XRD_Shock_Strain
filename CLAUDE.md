# xrd_profile — project-local conventions

This file captures conventions specific to the `xrd_profile` package
that future sessions in this directory should pick up automatically.
It does not duplicate the global scientific-writing / manuscript
style guide at `~/.claude/CLAUDE.md` — that still applies on top of
whatever is here.

## What this is

A Python toolkit for quantitative analysis of powder X-ray diffraction
peak profiles (Williamson-Hall, Warren-Averbach, Scherrer, PDF) with
reference-guided peak detection. Methodological subject of a Journal
of Applied Crystallography (JAC) manuscript and the analysis backbone
for a follow-on Meteoritics & Planetary Science (MAPS) study of HED
shock systematics.

Currently shipped: v0.3.0 (Phase 1, 2026-05-05). Next planned tags:
v0.4.0 → v0.5.0 → v1.0.0 (Phase 2). See spec / plan locations below.

## Working directory, Python environment, and test command

- Package root: `C:\Users\Matthew Izawa\Documents\Dan Applin\Llunr\xrd_profile\`.
- The package root is itself a git repo. The parent directory `Llunr/`
  is not a git repo; it holds the JAC manuscript, analysis pipeline
  scripts, and instrument data alongside the package.
- **Python interpreter**: Anaconda base env at
  `C:\Users\Matthew Izawa\anaconda3\python.exe` (Python 3.13.9). The
  `python` / `pip` / `pytest` commands are *not* on PATH from
  bash / PowerShell unless an Anaconda shell-init was run; invoke the
  interpreter by full path. From Git Bash:
  `"/c/Users/Matthew Izawa/anaconda3/python.exe"`.
  From PowerShell: `& "C:\Users\Matthew Izawa\anaconda3\python.exe"`.
- **Test command** (from the package root or any worktree):
  ```
  "/c/Users/Matthew Izawa/anaconda3/python.exe" -m pytest tests/ -v
  ```
  Stop on first failure with `-x`. The local `xrd_profile/` package
  becomes importable because pytest's rootdir-on-sys.path behaviour
  picks up the source tree — no `pip install -e .` needed.
- The base env already has numpy, scipy, matplotlib, pymatgen, pytest
  (verified 2026-05-05 against v0.3.0; the full v0.3.0 test suite
  passes at 59/59 in ~30 s).

## Branching, worktrees, and tags

- **Feature work happens in `.worktrees/<branch>/`** (project-local,
  hidden, gitignored). Branch naming: `feature/<release-version>` for
  release-bound work (e.g., `feature/v0.4.0`); `fix/<short-name>` for
  patches. Set up via `superpowers:using-git-worktrees`.
- **`main` holds tagged releases only.** Tags follow semver:
  `v0.2.0`, `v0.3.0`, `v0.4.0`, ... Each tag created on the
  release-complete commit on `main` after merging the feature branch.
- **Never push to remote without explicit user approval.** The remote
  is `https://github.com/MatthewIzawa/XRD_Shock_Strain` (currently
  private; public release deferred behind the JAC manuscript decision).
- **Never `git push --force` or `git rebase` published commits.** The
  v0.2.0 / v0.3.0 tags and their commits are reproducibility anchors
  for the JAC manuscript and the v0.3.0 release.

## Strict-additive contract and tiered goldens

- Every release tag (v0.2.0, v0.3.0, v0.4.0, ...) ships a
  `tests/fixtures/golden_v<version>_results.json` snapshot. Later
  releases retain numerical equality on every key already present in
  earlier-tier goldens (key-subset value-equality semantics —
  additional new keys at later tags are allowed, existing key values
  are not).
- `tests/test_backward_compat.py` exercises every active tier each
  release. Adding a tier never modifies an earlier tier's golden file.
- Regenerating an existing golden requires explicit reasoning in the
  commit message — counts as a numerical break and is rare. Use
  `python scripts/regenerate_goldens.py --tier v0.X.0`.

## Optional dependency extras

- `pip install xrd_profile` — core only (numpy, scipy, matplotlib).
- `pip install xrd_profile[cif]` — adds pymatgen for the `Phase` API.
- `pip install xrd_profile[cli]` — *(v1.0.0)* adds pyyaml + pydantic
  for the YAML CLI frontend.
- `pip install xrd_profile[docs]` — *(v1.0.0)* adds mkdocs-material +
  mkdocstrings for documentation builds.
- `pip install -e .[dev,cif,cli,docs]` — full contributor install.

The core install must stay minimal: every Stokes-deconvolution / W-H /
W-A / Scherrer / PDF code path runs on numpy arrays only.

## JAC pipeline isolation

The following paths under `Llunr/` are **off-limits** to package work
in any session — they are the JAC manuscript and its analysis pipeline,
which must not drift while the manuscript is in author / reviewer
review:

- `Llunr/Paper1_JAC/` (manuscript draft, figures, supplementary)
- `Llunr/paper1_figures.py`
- `Llunr/run_guided_wh.py`
- `Llunr/compile_survey.py`

If a package change requires a corresponding pipeline-script update,
flag it explicitly to the user — do not edit the pipeline scripts
unsolicited.

Also: `Llunr/xrd_profile/examples/legacy/` (verbatim v0.2.0 originals
preserved as historical artifacts) is similarly frozen.

## Spec / plan / changelog locations

- Design specs: `docs/superpowers/specs/YYYY-MM-DD-<topic>-design.md`.
  Active specs:
    - `2026-05-05-xrd-profile-v1-phase1-design.md` (v0.3.0, shipped).
    - `2026-05-05-xrd-profile-v1-phase2-design.md` (v0.4 / v0.5 / v1.0,
       approved, in progress).
- Implementation plans: `docs/superpowers/plans/YYYY-MM-DD-<feature>.md`.
  Active plans:
    - `2026-05-05-xrd-profile-v0.3.0-phase1.md` (executed).
    - `2026-05-05-xrd-profile-v0.4.0.md` (next to execute).
- Changelog: `CHANGELOG.md`. Every release tag has an entry in
  Keep-a-Changelog format.

## Commit-message style

Short imperative subject lines, ~50–60 characters, no leading verb-tag
(no `feat:` / `fix:` prefixes). Examples from existing log:
- `Release v0.3.0`
- `Update README to v0.3.0: Phase API, run_all, [cif] install`
- `Add [cif] optional extra for pymatgen, bump pyproject version to 0.3.0`

Co-author trailer is added by the harness — do not duplicate manually.

## When the user asks to ship a release

1. Confirm all plan tasks complete with passing tests.
2. Confirm `CHANGELOG.md` has an entry for the new tag.
3. Confirm `__version__` in `xrd_profile/__init__.py` and `pyproject.toml`
   match the new tag.
4. Merge the feature branch into `main` (fast-forward when possible).
5. Tag on `main` with `git tag -a vX.Y.Z -m "Release vX.Y.Z: <one-line summary>"`.
6. **Do not push** the tag to the remote without explicit user approval.

## Quick orientation pointers for new sessions

- Latest release: see `git tag --list 'v*' --sort=-v:refname | head -1`.
- What landed in the latest release: top entry of `CHANGELOG.md`.
- What's planned next: most recent plan in `docs/superpowers/plans/`.
- What architectural commitments are in force: most recent spec in
  `docs/superpowers/specs/`.
