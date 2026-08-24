# RaMAx documentation

This directory contains the detailed documentation for RaMAx 1.0.7.

## Getting started

- [Usage guide](usage.md): input format, MAF, HAL, and PAF commands, optimization
  behavior, and complete examples.
- [Command-line parameters](parameters.md): authoritative public CLI reference.

## Operation

- [Restart and work directories](restart.md): restart schema, retained state,
  and compatibility rules.
- [Output formats](output-formats.md): MAF/HAL/PAF/GFA conventions, validation, and
  normalized comparisons.
- [PAF-to-GFA workflow](paf-to-gfa-workflow.md): companion FASTA generation,
  direct seqwish induction, PGGB post-processing, validation, and a completed
  Chr09 example.
- [Troubleshooting](troubleshooting.md): common startup, dependency, resource,
  restart, and output problems.

## Internals

- [Realignment architecture](realignment_architecture.md): execution order,
  graph transactions, external MSA handling, and export consistency.
- [Changelog](CHANGELOG.md): release-level changes.

The root [README](../README.md) remains the short installation and quick-start
entry point. Run `ramax --help` to inspect the options in the installed binary.
