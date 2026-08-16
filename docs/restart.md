# Restart and work directories

## Configuration schema

RaMAx 1.0.5 writes all restartable user settings to:

```text
<work>/config.json
```

The root object is `common_args` and contains `schema_version: 1`. It records
the input/output paths, thread count, anchor and graph parameters, output
format, optimization flags and spans, logging settings, and HAL root/reference
selection.

PAF pair-selection mode is stored separately in `<work>/paf_mode.txt` so the
schema remains compatible with existing MAF/HAL work directories. Restart
uses the saved value; an older PAF work directory without this sidecar defaults
to `connected`.

## Starting a new run

Use a new or empty work directory:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -o /data/project/results/alignment.maf \
  -w /data/project/work/ramax-run \
  -t 16
```

In a Release build, an existing non-empty work directory is rejected. This
prevents accidental reuse of indexes or graph state from a different run.

RaMAx removes the entire work directory after a successful export. Keep the
final output outside it.

## Restarting an interrupted run

If a run fails or is interrupted after creating its work state, restart with:

```bash
ramax --restart -w /data/project/work/ramax-run
```

RaMAx loads the saved paths and parameters from `config.json`. Do not repeat
`--input`, `--output`, `--threads`, optimization options, or other algorithm
settings on the restart command. The saved configuration is authoritative.

Restart reuses compatible indexes and intermediate files; it does not promise
to resume inside an individual alignment operation that was interrupted.

## Compatibility

Only schema version 1 is accepted by RaMAx 1.0.5. Older experimental work
directories that lack `schema_version` or use a different schema are rejected
with an explicit incompatibility error. They are not silently upgraded or
interpreted with current defaults.

To rerun an old job, create a fresh work directory and provide the original
seqfile and desired current parameters as a new command.

## Files useful after failure

Before deleting a failed work directory, retain:

- `config.json` for the exact effective settings;
- `paf_mode.txt` for the effective PAF pair-selection mode, when present;
- `RaMAx.log` for the last completed stage and error;
- index and graph state when attempting `--restart`;
- `minipoa_tmp/` only when diagnosing a forced termination, because normal
  minipoa invocations clean their own temporary files.
