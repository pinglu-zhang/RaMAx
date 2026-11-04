#!/usr/bin/env python3
"""Check release metadata without executing build or upload hooks."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def read_contract(root, tag=''):
    cmake = (root / 'CMakeLists.txt').read_text(encoding='utf-8')
    recipe = (root / 'recipe/meta.yaml').read_text(encoding='utf-8')
    docker = (root / 'Dockerfile').read_text(encoding='utf-8')
    patterns = [(cmake, r'project\(RaMAx VERSION ([0-9]+\.[0-9]+\.[0-9]+)'),
                (recipe, r'{% set version = "([0-9]+\.[0-9]+\.[0-9]+)" %}'),
                (docker, r'(?m)^ARG RAMAX_VERSION=([0-9]+\.[0-9]+\.[0-9]+)$')]
    versions = []
    for text, pattern in patterns:
        matches = re.findall(pattern, text)
        if len(matches) != 1:
            raise ValueError('missing or ambiguous release version')
        versions.append(matches[0])
    if len(set(versions)) != 1:
        raise ValueError('CMake, Conda and Docker versions disagree')
    version = versions[0]
    if tag and tag != 'v' + version:
        raise ValueError('release tag does not match the source version')
    if 'halAppendCactusSubtree' in recipe or 'halAppendCactusSubtree' in docker:
        raise ValueError('obsolete Cactus append dependency in native HAL package')
    return version


def verify_source_tree(root):
    baseline = json.loads((root / 'ci/source-baseline.json').read_text(encoding='utf-8'))
    listing = subprocess.check_output(['git', 'ls-tree', '-d', 'HEAD', '--',
                                      *baseline['trees']], cwd=root, text=True)
    actual = {line.split('\t')[1]: line.split()[2] for line in listing.splitlines()}
    if actual != baseline['trees']:
        raise ValueError('production source differs from the frozen dev baseline')
    return baseline['commit']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, default=Path('.'))
    p.add_argument('--tag', default='')
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    version = read_contract(a.root, a.tag)
    baseline = verify_source_tree(a.root)
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=a.root, text=True).strip()
    result = {'version': version, 'commit': commit, 'baseline': baseline, 'tag': a.tag,
              'platform': 'linux-64 / linux/amd64', 'status': 'PASS'}
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(result, indent=2) + '\n')
    if os.environ.get('GITHUB_OUTPUT'):
        with open(os.environ['GITHUB_OUTPUT'], 'a') as f:
            f.write(f'version={version}\ncommit={commit}\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
