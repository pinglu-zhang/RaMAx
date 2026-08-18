# Restart and work directories

## What restart means

RaMAx 1.0.7 uses `--restart` to reuse expensive input preparation and
FM-index caches. It is not an alignment checkpoint system.

Reusable state:

- copied or downloaded raw FASTA;
- normalized alignment FASTA;
- HAL softmask sidecars;
- FM-index main, sampled-SA, and wavelet-tree files.

Anchor search, clustering, extension, DP, Block construction, graph
optimization, and every requested export always run again from the beginning.
RaMAx does not serialize these large in-memory objects.

## Starting and restarting

A new Release run requires a new or empty work directory:

```bash
ramax \
  -i /data/project/seqfile.txt \
  -w /data/project/work/ramax-run \
  -o /data/project/results/alignment.maf \
  -o /data/project/results/alignment.paf \
  -t 16
```

After interruption or export failure, reuse the work directory:

```bash
ramax --restart -w /data/project/work/ramax-run
```

Saved values are the defaults for restart. Explicit options override them:

```bash
ramax --restart -w /data/project/work/ramax-run \
  -t 24 \
  --min_anchor_length 30 \
  -o /data/project/results/retry.maf \
  -o /data/project/results/retry.paf \
  --paf-mode connected
```

`-i` cannot be used with `--restart`. If at least one `-o` is present, the
complete new repeated-output list replaces the saved list. Otherwise all saved
outputs are retained. Existing boolean flags keep their current one-way
behavior; for example `--one-round` can enable one-round mode but there is no
new restart-only flag to disable it.

The merged configuration is validated against the final output set. An
explicit `--paf-mode` requires PAF, an explicit `--gfa-version` requires GFA,
and an explicit `--root` requires HAL.
After input validation, the latest effective configuration is saved atomically
and becomes the default for a later restart.

## Configuration and input identity

`<work>/config.json` uses `schema_version: 4` and contains the complete output
list and all public alignment, graph, export, performance, and logging values.
`<work>/input_manifest.json` records the seqfile, species mappings, and the
size/mtime of local inputs.

The seqfile path and contents, species mappings, and local FASTA metadata must
remain unchanged. A mismatch stops before preprocessing and instructs the user
to start with a new work directory. A valid cached URL snapshot can be reused
without contacting the URL; the URL is accessed again only when that snapshot
is absent or invalid.

RaMAx reads schema-1, schema-2, and schema-3 work directories. Parameters absent
from older schemas use their compatibility defaults; in particular GFA version
migrates to `1.1`. Legacy schema-1 `outputs.json`,
`paf_mode.txt`, and `accurate_skip_threshold.txt` are read when present.
Existing legacy caches are trusted once, completion metadata is generated, and
the effective configuration is migrated to schema 4 with a warning.

## Cache validation and lifecycle

New raw, clean, softmask, and FM-index artifacts are written through partial
files and published only after successful completion. Their completion markers
record cache format, source identity, file size, and mtime. A missing,
incomplete, changed, or unloadable cache item is rebuilt independently.

FM indexes are still created lazily at the beginning of each reference round.
An index already completed by an earlier attempt is reused; missing later-round
indexes are built when reached. `--slow-build` affects only indexes that need to
be built or rebuilt.

Before restarted alignment begins, RaMAx:

1. archives the previous log as `RaMAx.restart.N.log`;
2. retains `data/`, `index/`, configuration, and completion markers;
3. removes stale `result/`, `mask_interval/`, and `minipoa_tmp/` directories;
4. logs that alignment is restarting from the beginning.

After every requested export succeeds, the complete work directory is removed.
If any exporter fails, other formats are still attempted, successful output
files are retained, the process exits nonzero, and the work directory remains
available for another restart.

## Useful failure evidence

Before manually deleting a failed work directory, retain:

- `config.json` and `input_manifest.json`;
- `RaMAx.log` and any `RaMAx.restart.N.log` archives;
- cache completion markers next to raw, clean, softmask, and FM-index files;
- the exact binary version and restart command.

`result/`, `mask_interval/`, and `minipoa_tmp/` are diagnostic artifacts only;
they are deliberately not used as restart checkpoints.
