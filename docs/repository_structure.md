# Repository structure

This repository is intended for sanitized code, configuration templates,
environment specifications, and documentation for a cross-cohort paediatric IBD
metagenomic study.

| Directory | Intended purpose |
| --- | --- |
| `scripts/server/` | Bash scripts for server-side metagenomic preprocessing, with Linux-compatible LF line endings and configurable paths. |
| `scripts/local/` | R scripts for local statistical analysis, machine-learning and random-forest modelling. |
| `config/` | Public configuration templates with clearly marked placeholders for machine-specific settings. |
| `environment/` | Software requirements and environment specifications for reproducibility, recorded after review. |
| `docs/` | Documentation of repository organization and reviewed workflows. |

Empty directories are not tracked by Git. They will become visible in clones
when reviewed files are added; placeholder files are not required.

Do not store raw or processed data, patient-level metadata, intermediate files,
model objects, figures, tables, analysis results, credentials, or private storage
information in these directories. Any specific public example requires explicit
user approval.
