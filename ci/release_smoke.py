#!/usr/bin/env python3
"""Small deterministic installed-product test. Never uses production inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import subprocess
import time


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def reverse(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def check_gfa(path, sequences):
    nodes, links, walks = {}, set(), []
    for line in path.read_text().splitlines():
        f = line.split('\t')
        if f[0] == 'S':
            require(f[1] not in nodes and f[2] != '*', 'duplicate/absent GFA sequence')
            nodes[f[1]] = f[2]
        elif f[0] == 'L':
            require(f[5] == '0M', 'unexpected GFA overlap')
            a, b = (f[1], f[2]), (f[3], f[4])
            links.add((a, b))
            links.add(((b[0], '+' if b[1] == '-' else '-'),
                       (a[0], '+' if a[1] == '-' else '-')))
        elif f[0] == 'W':
            walks.append(f)
    require(nodes and walks, 'GFA has no nodes or walks')
    require(all(a[0] in nodes and b[0] in nodes for a, b in links), 'dangling GFA link')
    recovered = {}
    for f in walks:
        steps = re.findall(r'([><])([^><]+)', f[6])
        require(''.join(d + n for d, n in steps) == f[6], 'malformed GFA walk')
        route = [(n, '+' if d == '>' else '-') for d, n in steps]
        require(all(n in nodes for n, _ in route), 'walk references absent node')
        require(all((a, b) in links for a, b in zip(route, route[1:])), 'missing walk edge')
        dna = ''.join(nodes[n] if d == '+' else reverse(nodes[n]) for n, d in route)
        key = (f[1], f[3])
        require(key not in recovered and int(f[4]) == 0 and int(f[5]) == len(dna),
                'duplicate or partial GFA walk')
        recovered[key] = dna
    require(recovered == sequences, 'GFA does not reconstruct every complete input contig')
    return {'nodes': len(nodes), 'links_bidirected': len(links), 'walks': len(walks)}


def check_maf(path, sequences):
    expected = {s + '.' + c: dna for (s, c), dna in sequences.items()}
    rows, matching_bases = 0, 0
    for line in path.read_text().splitlines():
        if not line.startswith('s '):
            continue
        f = line.split()
        require(len(f) == 7, 'malformed MAF row')
        if f[1] not in expected:  # Internal ancestors may also be represented.
            continue
        dna = expected[f[1]]
        start, size, length = int(f[2]), int(f[3]), int(f[5])
        require(length == len(dna) and 0 <= start <= start + size <= length, 'MAF coordinates')
        oriented = dna if f[4] == '+' else reverse(dna)
        require(f[4] in ('+', '-') and f[6].replace('-', '').upper() ==
                oriented[start:start + size].upper(), 'MAF DNA mismatch')
        rows += 1
        matching_bases += size
    require(rows > 2 and matching_bases > 30000, 'MAF contains no substantial alignment')
    return {'leaf_rows': rows, 'leaf_row_bases': matching_bases}


def check_paf(path, sequences):
    lengths = {s + '.' + c: len(dna) for (s, c), dna in sequences.items()}
    count = 0
    for line in path.read_text().splitlines():
        f = line.split('\t')
        require(len(f) >= 12 and f[0] in lengths and f[5] in lengths, 'unknown PAF identity')
        qlen, qs, qe, tlen, ts, te = map(int, (f[1], f[2], f[3], f[6], f[7], f[8]))
        require(qlen == lengths[f[0]] and tlen == lengths[f[5]], 'PAF sequence lengths')
        require(0 <= qs < qe <= qlen and 0 <= ts < te <= tlen and f[4] in ('+', '-'),
                'PAF coordinates')
        tags = [x[5:] for x in f[12:] if x.startswith('cg:Z:')]
        require(len(tags) == 1, 'PAF lacks a unique CIGAR')
        parts = re.findall(r'(\d+)([MID=X])', tags[0])
        require(''.join(n + c for n, c in parts) == tags[0], 'invalid PAF CIGAR')
        require(sum(int(n) for n, c in parts if c != 'D') == qe - qs and
                sum(int(n) for n, c in parts if c != 'I') == te - ts, 'PAF CIGAR consumption')
        count += 1
    require(count > 0, 'empty PAF')
    return {'records': count}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--ramax', default='ramax')
    p.add_argument('--docker-image')
    p.add_argument('--hal-validator', type=Path, required=True)
    p.add_argument('--version', default='1.1.0')
    args = p.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    result = {'status': 'RUNNING', 'version': args.version, 'commands': []}
    def save():
        (root / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
    save()
    try:
        rng = random.Random(11002026)
        dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        ref, aux = dna(60000), dna(3001)
        a = list(ref)
        for i in range(903, len(a), 997):
            a[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[a[i]]
        a = ''.join(a)
        a = a[:18000] + dna(137) + a[18000:41000] + a[41113:]
        b = ref[:22000] + ref[22081:47000] + dna(171) + ref[47000:]
        sequences = {(sample, contig): value for sample, primary in [('R', ref), ('A', a), ('B', b)]
                     for contig, value in [('chr1', primary), ('aux', aux)]}
        seqfile = root / 'seqfile.tsv'
        entries = ['((R:1,A:1)ancRA:1,B:1)root;\n']
        for sample in ('R', 'A', 'B'):
            fasta = root / (sample + '.fa')
            fasta.write_text(''.join('>' + c + '\n' + sequences[(sample, c)] + '\n'
                                     for c in ('chr1', 'aux')))
            entries.append(sample + '\t' + str(fasta) + '\n')
        seqfile.write_text(''.join(entries))
        expected = root / 'expected.tsv'
        expected.write_text(''.join(s + '\t' + c + '\t' + v + '\n'
                                    for (s, c), v in sequences.items()))
        if args.docker_image:
            prefix = ['docker', 'run', '--rm', '--user', f'{os.getuid()}:{os.getgid()}',
                      '-v', f'{root}:{root}', '-w', str(root), args.docker_image]
        else:
            prefix = [args.ramax]
        version = subprocess.check_output(prefix + ['--version'], text=True, timeout=60).strip()
        require(version == 'RaMAx version ' + args.version, 'installed version mismatch')
        command = prefix + ['-i', str(seqfile), '-w', str(root / 'work'), '-t', '2', '--ref', 'R']
        for suffix in ('maf', 'hal', 'paf', 'gfa'):
            command += ['-o', str(root / ('alignment.' + suffix))]
        result['commands'].append(command)
        save()
        env = os.environ.copy()
        env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')
        started = time.monotonic()
        with (root / 'run.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env,
                           timeout=600, check=True)
        result['wall_seconds'] = time.monotonic() - started
        result['gfa'] = check_gfa(root / 'alignment.gfa', sequences)
        result['maf'] = check_maf(root / 'alignment.maf', sequences)
        result['paf'] = check_paf(root / 'alignment.paf', sequences)
        validator = [str(args.hal_validator.resolve()), str(root / 'alignment.hal'), str(expected)]
        result['commands'].append(validator)
        result['hal'] = subprocess.check_output(validator, text=True, timeout=120).strip()
        result['input_sha256'] = {s + '/' + c: hashlib.sha256(v.encode()).hexdigest()
                                  for (s, c), v in sequences.items()}
        result['outputs'] = {x.name: {'bytes': x.stat().st_size,
                                     'sha256': hashlib.sha256(x.read_bytes()).hexdigest()}
                             for x in root.glob('alignment.*')}
        result['status'] = 'PASS'
        save()
    except BaseException as e:
        result.update(status='FAILED', error=repr(e))
        save()
        raise


if __name__ == '__main__':
    main()
