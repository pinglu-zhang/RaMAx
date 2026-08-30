# Restart and work directories

## What restart means

RaMAx 1.0.8 and later use `--restart` to reuse expensive input preparation. The
suffix-array index is memory-only and is rebuilt by every process; restart is
not an alignment checkpoint system.

Reusable state:

- copied or downloaded raw FASTA;
- normalized alignment FASTA;
- HAL softmask sidecars;

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
  --sa-sampling-rate 1 \
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
an explicit `--gfa-profile` requires GFA, and an explicit `--root` requires HAL.
After input validation, the latest effective configuration is saved atomically
and becomes the default for a later restart.

## Configuration and input identity

`<work>/config.json` uses `schema_version: 6` and contains the complete output
list and all public alignment, graph, export, performance, and logging values.
`<work>/input_manifest.json` records the seqfile, species mappings, and the
size/mtime of local inputs.

The seqfile path and contents, species mappings, and local FASTA metadata must
remain unchanged. A mismatch stops before preprocessing and instructs the user
to start with a new work directory. A valid cached URL snapshot can be reused
without contacting the URL; the URL is accessed again only when that snapshot
is absent or invalid.

RaMAx reads schema-1 through schema-5 work directories. Parameters absent from
older schemas use their compatibility defaults; in particular schema 1-5
workdirs migrate to `--sa-sampling-rate 1`, GFA profile migrates to `exact`,
and schemas before 4 migrate GFA version to `1.1`. Legacy schema-1 `outputs.json`,
`paf_mode.txt`, and `accurate_skip_threshold.txt` are read when present.
Existing legacy caches are trusted once, completion metadata is generated, and
the effective configuration is migrated to schema 6 with a warning.

## Cache validation and lifecycle

New raw, clean, and softmask artifacts are written through partial files and
published only after successful completion. Their completion markers record
cache format, source identity, file size, and mtime. A missing, incomplete,
changed, or unloadable cache item is rebuilt independently.

Suffix-array indexes are created in memory at the beginning of each reference round and are never read from or written to disk. The persisted `--sa-sampling-rate` value must be `1`; other values are rejected before index construction. Restart therefore rebuilds the suffix array while continuing to reuse eligible preprocessing artifacts. Existing `.saidx` files from older runs are ignored and are not deleted.
Each reference round builds its suffix array when reached. `--slow-build`
continues to control the historical text/sentinel layout; it does not enable a
disk cache.

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
- cache completion markers next to raw, clean, and softmask files;
- the exact binary version and restart command.

`result/`, `mask_interval/`, and `minipoa_tmp/` are diagnostic artifacts only;
they are deliberately not used as restart checkpoints.
