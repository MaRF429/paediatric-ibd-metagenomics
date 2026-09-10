# Repository guidance

## Repository scope
- This repository contains reproducible code, configuration templates, environment specifications, and documentation for a cross-cohort paediatric IBD metagenomic study.
- Do not add raw data, processed data, patient-level metadata, intermediate files, model objects, figures, tables, or analysis results unless the user explicitly approves a specific public example.
- Public dataset accession numbers may be included, but personal identifiers, credentials and private storage information must never be included.

## Source-file protection
- Files outside this Git repository are source materials and must be treated as read-only.
- Never edit, rename, move or delete original server scripts or original local R scripts.
- When source scripts are provided, create sanitized copies inside this repository.

## Scientific integrity
- Preserve analytical parameters, statistical methods, random seeds, feature definitions, cohort assignments, outcomes and comparison directions.
- Do not silently change scientific logic to make a script cleaner.
- Explicit user approval is required before changing algorithms, thresholds, covariates, feature-selection rules, validation design or model-training logic.
- Discovery cohorts may be used for feature selection, tuning and model development.
- Independent validation cohorts must only be used for locked-model evaluation.
- External validation results must not influence feature selection, algorithm ranking or hyperparameter tuning.
- Preserve the five-algorithm feature-importance workflow, the 113 algorithmic combinations and the seven-stage random-forest workflow unless the user explicitly requests a change.

## Privacy and portability
- Remove absolute local paths, server usernames, account names, home-directory paths, storage mount points, API keys, tokens, passwords and private URLs from repository copies.
- Replace machine-specific paths with command-line arguments, configuration files or clearly marked placeholders.
- Do not expose sample identifiers or individual-level clinical metadata.
- Never print or inspect credentials unnecessarily.

## Shell scripts
- Server scripts must use Bash and Linux-compatible LF line endings.
- Prefer `set -euo pipefail` when compatible with the original workflow.
- Quote path variables and validate required inputs.
- Preserve software parameters and database-version settings.
- Do not execute computationally intensive metagenomic pipelines on the local computer.
- Validate shell syntax with `bash -n` when Bash is available.

## R scripts
- Repository-facing R scripts and comments should be written in English.
- Use project-relative paths or configuration templates rather than absolute paths.
- Do not automatically install packages inside analysis scripts.
- Record package requirements separately.
- Preserve random seeds and data-splitting boundaries.
- Use parse-only validation when real input data are unavailable.
- Do not execute full analyses unless the user explicitly requests it.

## Git safety
- Work on task-specific branches, not directly on `main`.
- Do not commit, push, merge, create releases or change repository visibility without explicit user approval.
- Do not use `git add .`; stage only specifically reviewed files.
- Do not use destructive Git commands.
- Before proposing a commit, report every modified or untracked file.

## Communication
- Explain work and warnings to the user in Chinese.
- Keep repository files, code comments and public documentation in English.
- If a scientific or privacy-sensitive decision is uncertain, stop and ask the user instead of guessing.
