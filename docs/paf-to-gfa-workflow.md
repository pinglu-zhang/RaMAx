# RaMAx PAF 到泛基因组 GFA

本文说明如何从一个 RaMAx `seqfile` 出发，生成名称严格对应的合并
FASTA 和 connected PAF，再通过 seqwish 或 PGGB 获得 GFAv1 泛基因组图。

适用版本：`dev-pan 1.4`。

## 工作流概览

有两条可用路线：

1. **seqwish 最短路线**：直接将 RaMAx PAF 诱导为原始 GFA，适合首先验证
   PAF、路径名称和图结构是否正确。
2. **PGGB 后处理路线**：通过 `pggb -a` 跳过 wfmash，复用 RaMAx PAF，
   再执行 seqwish、smoothxg、gfaffix 和 odgi 处理，适合正式分析。

建议先跑通 seqwish 路线，再运行 PGGB。两条路线都必须使用由同一个
seqfile 生成的配套 FASTA，不能手工改写 FASTA 或 PAF 中的名称。

## 前置条件

需要以下命令：

- `ramax`
- `ramax-paf-fasta`
- seqwish 路线：`seqwish`，推荐同时安装 `odgi` 或 `vg`
- PGGB 路线：`pggb`、`samtools`，以及 PGGB 自带的下游依赖

检查版本：

```bash
ramax --version
ramax-paf-fasta --version
seqwish --version
pggb --version
odgi version
vg version
```

`ramax` 与 `ramax-paf-fasta` 应来自同一份 RaMAx 源码或同一个安装前缀。

## 输入名称约定

PAF 中的序列名称固定为：

```text
species.original_contig_name
```

`original_contig_name` 是原 FASTA header 的第一个空白字段。工具不会自动
删除已经存在的物种前缀，也不会为重复名称增加后缀。

例如 seqfile 为：

```text
ref       /data/ref.fa
sample1   /data/sample1.fa.gz
sample2   /data/sample2.fa
```

若三个 FASTA 都包含 `>Chr09`，输出名称为：

```text
ref.Chr09
sample1.Chr09
sample2.Chr09
```

PAF 不要求 Newick 树；seqfile 带有可选 Newick 首行时，
`ramax-paf-fasta` 也可以识别。

## 建议的实验目录

代码版本、实验版本和结果目录应分开记录。下面的命令使用独立的 1.4
目录，并把临时文件放在数据盘：

```bash
set -euo pipefail

RAMAX=/path/to/ramax
RAMAX_PAF_FASTA=/path/to/ramax-paf-fasta
SEQFILE=/data/project/genomes.seqfile
RUN=/mnt/d/Result/RaMAx/exp/pangenome/my_dataset_devpan_1.4
THREADS=12

mkdir -p "$RUN/input" "$RUN/output" "$RUN/logs" "$RUN/tmp/seqwish"

"$RAMAX" --version | tee "$RUN/logs/ramax.version.txt"
"$RAMAX_PAF_FASTA" --version \
  | tee "$RUN/logs/ramax-paf-fasta.version.txt"
seqwish --version | tee "$RUN/logs/seqwish.version.txt"
```

RaMAx 新运行的 workdir 必须为空。运行成功后，RaMAx 可能自动删除该
workdir。最终 PAF、FASTA 和日志不要放进 RaMAx workdir。

## 路线 A：直接使用 seqwish

### 1. 生成配套 FASTA

seqwish 支持普通或 gzip FASTA，因此可以直接生成 `.fa.gz`：

```bash
"$RAMAX_PAF_FASTA" \
  -i "$SEQFILE" \
  -o "$RUN/input/sequences.fa.gz"
```

该工具会：

- 按 seqfile 中的物种顺序和原 FASTA contig 顺序输出；
- 将 header 改为 `species.contig`；
- 将序列转为大写；
- 保留 `A/T/G/C`，将其他字符转为 `N`；
- 在全部检查通过后原子发布输出文件。

### 2. 生成 connected PAF

```bash
"$RAMAX" \
  -i "$SEQFILE" \
  -o "$RUN/output/alignment.connected.paf" \
  -w "$RUN/work-ramax" \
  -t "$THREADS" \
  --paf-mode connected \
  --log-level info \
  2>&1 | tee "$RUN/logs/ramax.log"
```

运行结束后检查：

```bash
grep 'PAF export complete' "$RUN/logs/ramax.log"
```

正式验收要求 `invalid=0`。`fallback_blocks` 最好为 0；如果不为 0，表示
个别 Block 已退化为 all2all 输出，仍可生成 PAF，但应在实验记录中说明。

### 3. 检查 FASTA 与 PAF 名称

```bash
gzip -dc "$RUN/input/sequences.fa.gz" \
  | awk '/^>/{sub(/^>/,""); split($0,a,/ /); print a[1]}' \
  | sort -u > "$RUN/logs/fasta.names.txt"

cut -f1,6 "$RUN/output/alignment.connected.paf" \
  | tr '\t' '\n' | sort -u > "$RUN/logs/paf.names.txt"

comm -3 "$RUN/logs/fasta.names.txt" "$RUN/logs/paf.names.txt"
```

如果 `comm` 输出只包含 FASTA 一侧的名称，可能是某条输入序列没有进入任何
有效 Block；应先确认是否符合实验预期。PAF 中出现 FASTA 不存在的名称则是
错误，不能继续运行 seqwish。

### 4. 生成原始 GFA

```bash
mkdir -p "$RUN/tmp/seqwish"

seqwish \
  -s "$RUN/input/sequences.fa.gz" \
  -p "$RUN/output/alignment.connected.paf" \
  -g "$RUN/output/graph.raw.gfa" \
  -b "$RUN/tmp/seqwish" \
  -t "$THREADS" \
  -k 0 \
  -P \
  2>&1 | tee "$RUN/logs/seqwish.log"
```

必须显式使用 `-k 0`。这里的 `-k` 是 seqwish 的最小精确匹配长度；大于
0 的值会主动删除短精确匹配，不再处于 RaMAx connected 模式的信息完整性
承诺范围内。

seqwish 使用磁盘支持的临时数据结构。务必用 `-b` 把临时目录放到空间充足
的数据盘，不要依赖当前目录或 `/home`。

### 5. 验证原始 GFA

```bash
test -s "$RUN/output/graph.raw.gfa"

odgi build \
  -g "$RUN/output/graph.raw.gfa" \
  -o "$RUN/output/graph.raw.og" \
  -s -O -t "$THREADS"

odgi validate \
  -i "$RUN/output/graph.raw.og" \
  -t "$THREADS"

odgi paths \
  -i "$RUN/output/graph.raw.og" \
  -L | sort > "$RUN/logs/graph.paths.txt"
```

若已安装 vg，也可以补充：

```bash
vg validate "$RUN/output/graph.raw.gfa"
```

到这里已经成功获得可用 GFA。smoothxg 和 gfaffix 属于图规范化步骤，不是
生成 GFA 的必要条件。

## 路线 B：接入 PGGB 后续流程

PGGB 的 `-a/--input-paf` 会跳过 wfmash，直接使用已有 PAF。之后 PGGB
继续运行 seqwish、smoothxg、gfaffix 和 odgi 相关步骤。

### 1. 生成普通 FASTA 并建立索引

PGGB 运行前要求 FASTA 存在 `.fai`。这里推荐输出未压缩 `.fa`：

```bash
"$RAMAX_PAF_FASTA" \
  -i "$SEQFILE" \
  -o "$RUN/input/sequences.fa"

samtools faidx "$RUN/input/sequences.fa"
```

`ramax-paf-fasta` 的 `.gz` 输出是普通 gzip，而 `samtools faidx` 对压缩
FASTA 通常要求 BGZF。直接使用普通 `.fa` 可以避免索引格式问题。

### 2. 准备 PAF

如果路线 A 已经使用相同 seqfile 生成了 PAF，可以直接复用：

```text
$RUN/output/alignment.connected.paf
```

不需要再运行 wfmash，也不要让 PGGB 重新替换 RaMAx 的同源关系。

### 3. 运行 PGGB

RaMAx 名称是 `species.contig`，不是 PanSN 格式，因此需要使用 `-n` 明确
输入物种或样本数。下面以 7 个输入为例：

```bash
N_GENOMES=7

mkdir -p "$RUN/pggb_out" "$RUN/tmp/pggb"

pggb \
  -i "$RUN/input/sequences.fa" \
  -o "$RUN/pggb_out" \
  -a "$RUN/output/alignment.connected.paf" \
  -n "$N_GENOMES" \
  -t "$THREADS" \
  -T 8 \
  -k 0 \
  -D "$RUN/tmp/pggb" \
  -v \
  2>&1 | tee "$RUN/logs/pggb.log"
```

参数含义：

- `-a`：使用 RaMAx PAF，跳过 wfmash；
- `-n`：输入物种或单倍型数量；
- `-t`：总线程数；
- `-T`：POA 线程数，内存紧张时应低于总线程数；
- `-k 0`：保持 connected PAF 的短精确匹配关系；
- `-D`：指定 PGGB 临时目录。

使用自定义 PAF 时，`-p` 和 `-s` 不再决定 wfmash 比对，因为 wfmash 已被
`-a` 跳过。初次跑通应尽量减少额外参数，再针对不同物种单独调整 PGGB 和
smoothxg 参数。

最终 GFA 通常为：

```bash
find "$RUN/pggb_out" \
  -maxdepth 1 \
  -name '*.final.gfa' \
  -type f \
  -print
```

### 4. 验证 final GFA

```bash
FINAL_GFA=$(find "$RUN/pggb_out" \
  -maxdepth 1 \
  -name '*.final.gfa' \
  -type f | head -n 1)

test -n "$FINAL_GFA"
test -s "$FINAL_GFA"

vg validate "$FINAL_GFA"

odgi build \
  -g "$FINAL_GFA" \
  -o "$RUN/output/graph.final.og" \
  -s -O -t "$THREADS"

odgi validate \
  -i "$RUN/output/graph.final.og" \
  -t "$THREADS"

odgi paths \
  -i "$RUN/output/graph.final.og" \
  -L | sort > "$RUN/logs/graph.final.paths.txt"
```

## 成功标准

一次完整运行应同时满足：

- `ramax-paf-fasta` 退出码为 0，qualified header 没有冲突；
- RaMAx 退出码为 0；
- `PAF export complete` 中 `invalid=0`；
- PAF 的 query/target 名称均存在于配套 FASTA；
- seqwish 或 PGGB 退出码为 0；
- 目标 GFA 非空；
- `odgi validate` 或 `vg validate` 通过；
- 图中的路径数和路径名与预期输入一致；
- 已保存命令、软件版本、日志以及输入、PAF 和 GFA 的 SHA-256。

## 常见问题

### seqwish 找不到序列名

FASTA header 与 PAF 名称不一致。必须用同一个 seqfile 运行
`ramax-paf-fasta`，不要手工重命名。

### PGGB 启动后立即提示缺少 `.fai`

输入 FASTA 没有索引。使用普通 `.fa` 输出，然后运行：

```bash
samtools faidx sequences.fa
```

### seqwish 临时文件占满当前目录

没有设置 `-b`，或临时目录所在分区空间不足。将 `-b` 指向数据盘上的独立
目录，并在确认任务失败后清理对应临时目录。

### connected 与 all 的下游图不一致

首先确认 seqwish 或 PGGB 使用了 `-k 0`，并检查 RaMAx PAF 汇总中的
`invalid` 和 `fallback_blocks`。大于 0 的 `-k` 会改变比对关系集合。

### PGGB 在 smoothxg 阶段内存不足

优先降低 `-T`，必要时再降低 `-t`。应根据日志确认失败发生在 smoothxg，
不要把所有 PGGB 失败都归因于 RaMAx PAF。

### RaMAx 提示 workdir 非空

新运行必须使用空 workdir。失败任务需要续跑时才使用 `--restart`；否则使用
新的实验目录，避免混合不同版本的中间状态。

## Chr09 真实流程验收记录

2026-08-17，使用 7 条高粱 Chr09 原始含 N 输入完成了
RaMAx PAF → PGGB final GFA 流程：

| 项目 | 结果 |
|---|---|
| 输入 | 7 条 qualified 路径，409,138,992 bp |
| RaMAx PAF | 16,366 个有效 Block，`invalid=0`，`fallback=0` |
| PAF 记录数 | 57,632 |
| PAF 大小 | 7,760,223 bytes |
| RaMAx 资源 | 4:09.45，最大 RSS 7,446,840 KB |
| PGGB | 0.7.4，使用 `-a` 跳过 wfmash，退出码 0 |
| PGGB 资源 | 30:13.04，最大 RSS 4,579,072 KB |
| final GFA | 104,772,632 bytes |
| 图结构 | 254,592 segments，345,752 links，7 paths |
| GFA 校验 | vg 1.74.0：`graph: valid` |
| GFA SHA-256 | `1da3ec92925c48405af3bcc4d8c7bf4588f27171637c9f93efffa3205833c67c` |

最终 GFA：

```text
/mnt/d/Result/RaMAx/exp/pangenome/Benchmarking_graph_pipelines_mini/results/
ramax_sorghum_sim_chr09_rawN_devpan_paf_pggb_20260817/pggb_out/
ramax_chr09.paf.2cec80e.c7b9bfa.smooth.final.gfa
```

该成功实例为了保持既有 PGGB benchmark 参数使用了 `pggb -k 25`。它证明
工程流程可以运行，但不属于 connected 模式严格的 `-k 0` 完整性验收。
正式 connected 实验应另跑 `-k 0`，并比较路径和碱基关系。

## 实验记录建议

```bash
sha256sum \
  "$SEQFILE" \
  "$RUN/input/sequences.fa" \
  "$RUN/output/alignment.connected.paf" \
  "$FINAL_GFA" \
  > "$RUN/logs/SHA256SUMS.txt"

grep 'PAF export complete' "$RUN/logs/ramax.log" \
  > "$RUN/logs/paf-export-summary.txt"
```

当前阶段统一标记为 `dev-pan 1.4`。若工作树包含未提交修改，除 commit 外还
应保存：

```bash
git status --short
git diff --stat
```

## 参考资料

- [seqwish usage](https://github.com/pangenome/seqwish#usage)
- [seqwish graph-induction algorithm](https://github.com/pangenome/seqwish#squish-graph-induction-algorithm)
- [PGGB](https://github.com/pangenome/pggb)
- [ODGI documentation](https://odgi.readthedocs.io/)
