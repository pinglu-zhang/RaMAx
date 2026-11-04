#!/usr/bin/env python3
"""Publish retained, validated artifacts. Never rebuild, retag Git or force uploads."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import urllib.error
import urllib.request

from check_release import read_contract, verify_source_tree

REPO = 'pinglu-zhang/RaMAx'
IMAGE = 'pingluzhang/ramax'


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def command(argv, capture=False):
    # In particular, never embed a token-bearing CalledProcessError in a log.
    r = subprocess.run(argv, check=False, text=True,
                       stdout=subprocess.PIPE if capture else None)
    if r.returncode:
        raise RuntimeError(f'{Path(argv[0]).name} failed with exit code {r.returncode}')
    return r.stdout if capture else ''


def gh_json(endpoint):
    return json.loads(command(['gh', 'api', endpoint], capture=True))


def unique(root, name):
    files = list(root.rglob(name))
    require(len(files) == 1, f'expected exactly one artifact {name}; found {len(files)}')
    return files[0]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def http_json(url, allow_missing=False):
    try:
        with urllib.request.urlopen(urllib.request.Request(
                url, headers={'User-Agent': 'RaMAx-release'}), timeout=60) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404 and allow_missing:
            return None
        raise RuntimeError(f'registry metadata lookup failed: HTTP {e.code}') from None


def remote_digest(url):
    h = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=120) as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--tag', required=True)
    p.add_argument('--artifacts', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--validated-run-id')
    a = p.parse_args()
    version = read_contract(Path('.'), a.tag)
    verify_source_tree(Path('.'))
    commit = command(['git', 'rev-parse', 'HEAD'], capture=True).strip()
    require(os.environ.get('GITHUB_REPOSITORY') == REPO, 'unexpected release repository')
    require(bool(os.environ.get('ANACONDA_API_TOKEN')), 'ANACONDA_API_TOKEN is required')
    a.work.mkdir(parents=True, exist_ok=False)
    state = {'tag': a.tag, 'version': version, 'commit': commit, 'status': 'RUNNING', 'stages': []}
    def stage(name):
        state['stages'].append(name)
        (a.work / 'STATUS.json').write_text(json.dumps(state, indent=2) + '\n')
    stage('VERIFYING_INPUTS')
    try:
        if a.validated_run_id:
            require(a.validated_run_id.isdecimal(), 'invalid validation run ID')
            run = gh_json(f'repos/{REPO}/actions/runs/{a.validated_run_id}')
            require(run['head_sha'] == commit and run['event'] == 'push' and
                    run['path'] == '.github/workflows/release.yml',
                    'resume artifacts must originate from the same tagged-source validation run')
        tag = gh_json(f'repos/{REPO}/git/ref/tags/{a.tag}')['object']
        if tag['type'] == 'tag':
            tag = gh_json(f'repos/{REPO}/git/tags/{tag["sha"]}')['object']
        require(tag['type'] == 'commit' and tag['sha'] == commit, 'Git tag/source mismatch')
        source = json.loads(unique(a.artifacts / 'source-evidence', 'source.json').read_text())
        require(source['commit'] == commit and source['version'] == version and
                source['status'] == 'PASS', 'source validation identity mismatch')
        for group, filename in [('source-evidence', 'source-smoke.json'),
                                ('conda-package', 'conda-smoke.json'),
                                ('docker-image', 'docker-smoke.json')]:
            data = json.loads(unique(a.artifacts / group, filename).read_text())
            require(data['status'] == 'PASS' and data['version'] == version,
                    f'{group} did not pass installed-product checks')
        package = unique(a.artifacts / 'conda-package', f'ramax-{version}-*.conda')
        package_sums = (package.parent.parent / 'SHA256SUMS').read_text().split()
        require(len(package_sums) == 2 and package_sums[1] == 'linux-64/' + package.name and
                package_sums[0] == digest(package), 'retained Conda checksum mismatch')
        tar = unique(a.artifacts / 'docker-image', 'image.tar.gz')
        sums = (tar.parent / 'SHA256SUMS').read_text().split()
        require(len(sums) == 2 and sums[1] == 'image.tar.gz' and sums[0] == digest(tar),
                'retained image checksum mismatch')
        expected_image = json.loads(unique(a.artifacts / 'docker-image', 'image-inspect.json').read_text())[0]
        require(expected_image['Config']['Labels']['org.opencontainers.image.revision'] == commit,
                'container source revision mismatch')
        require(expected_image['Config']['Labels']['org.opencontainers.image.version'] == version,
                'container version mismatch')
        command(['docker', 'load', '-i', str(tar)])
        image_id = expected_image['Id']
        loaded = json.loads(command(['docker', 'image', 'inspect', image_id], capture=True))[0]
        require(loaded['Id'] == image_id, 'loaded container identity mismatch')

        # Create or reuse a draft; an existing published release is never clobbered.
        releases = gh_json(f'repos/{REPO}/releases?per_page=100')
        for item in releases:
            match = re.fullmatch(r'v(\d+)\.(\d+)\.(\d+)', item['tag_name'])
            if match and not item['draft'] and not item['prerelease']:
                require(tuple(map(int, match.groups())) <= tuple(map(int, version.split('.'))),
                        'a newer stable release exists; refusing to change latest aliases')
        existing = next((x for x in releases if x['tag_name'] == a.tag), None)
        if existing is None:
            command(['gh', 'release', 'create', a.tag, '--repo', REPO, '--verify-tag',
                     '--draft', '--title', 'RaMAx ' + version,
                     '--notes-file', f'docs/releases/{a.tag}.md'])
        stage('DRAFT_READY')

        package_sha = digest(package)
        release = http_json(f'https://api.anaconda.org/release/malab/ramax/{version}', True)
        published = next((f for f in release.get('distributions', release.get('files', []))
                          if f['basename'] == 'linux-64/' + package.name), None) if release else None
        if published:
            url = published.get('download_url',
                f'https://api.anaconda.org/download/malab/ramax/{version}/linux-64/{package.name}')
            if url.startswith('//'):
                url = 'https:' + url
            require(remote_digest(url) == package_sha, 'different Conda artifact already published')
        else:
            command(['anaconda', '-t', os.environ['ANACONDA_API_TOKEN'], 'upload',
                     '--user', 'malab', str(package)])
        stage('CONDA_UPLOADED')

        version_image = IMAGE + ':' + version
        published_image = http_json(f'https://hub.docker.com/v2/repositories/{IMAGE}/tags/{version}/', True)
        if published_image:
            command(['docker', 'pull', version_image])
            remote_id = command(['docker', 'image', 'inspect', '--format', '{{.Id}}', version_image], True).strip()
            require(remote_id == image_id, 'different immutable container tag already published')
        else:
            command(['docker', 'tag', image_id, version_image])
            command(['docker', 'push', version_image])
        stage('DOCKER_VERSION_UPLOADED')

        env = a.work / 'public-conda'
        command(['conda', 'create', '-y', '-p', str(env), '--override-channels',
                 '-c', 'malab', '-c', 'bioconda', '-c', 'conda-forge',
                 'ramax=' + version, 'python=3.11'])
        validator = unique(a.artifacts / 'source-evidence', 'halValidate').resolve()
        validator.chmod(0o755)
        common = ['ci/release_smoke.py', '--version', version, '--hal-validator', str(validator)]
        command(['conda', 'run', '-p', str(env), 'python', *common,
                 '--ramax', str(env / 'bin/ramax'), '--output', str(a.work / 'public-conda-smoke')])
        command(['docker', 'pull', version_image])
        command(['python', *common, '--docker-image', version_image,
                 '--output', str(a.work / 'public-docker-smoke')])
        final_image = json.loads(command(['docker', 'image', 'inspect', version_image], True))[0]
        require(final_image['Id'] == image_id, 'public container identity changed')
        stage('PUBLIC_INSTALLATIONS_PASS')

        assets = a.work / 'assets'
        assets.mkdir()
        shutil.copyfile(package, assets / package.name)
        info = {'version': version, 'commit': commit, 'tag': a.tag,
                'conda': {'channel': 'malab', 'file': package.name, 'sha256': package_sha,
                          'platform': 'linux-64'},
                'docker': {'image': version_image, 'image_id': image_id,
                           'repo_digests': final_image['RepoDigests'], 'platform': 'linux/amd64'},
                'validation': {'source': 'PASS', 'conda': 'PASS', 'docker': 'PASS',
                               'public_installations': 'PASS'}}
        (assets / 'BUILD_INFO.json').write_text(json.dumps(info, indent=2, sort_keys=True) + '\n')
        (assets / 'SHA256SUMS').write_text(''.join(digest(f) + '  ' + f.name + '\n'
                                                for f in sorted(assets.iterdir())))
        remote_release = gh_json(f'repos/{REPO}/releases/tags/{a.tag}')
        for f in sorted(assets.iterdir()):
            old = next((x for x in remote_release['assets'] if x['name'] == f.name), None)
            if old:
                saved = a.work / ('existing-' + f.name)
                # gh streams binary assets directly; never decode a Conda file as text.
                with saved.open('wb') as output:
                    done = subprocess.run(['gh', 'api', f'repos/{REPO}/releases/assets/{old["id"]}',
                                           '-H', 'Accept: application/octet-stream'], stdout=output)
                require(done.returncode == 0 and digest(saved) == digest(f),
                        'existing Release asset differs; refusing overwrite')
            else:
                command(['gh', 'release', 'upload', a.tag, str(f), '--repo', REPO])
        for alias in ('.'.join(version.split('.')[:2]), 'latest'):
            command(['docker', 'tag', image_id, IMAGE + ':' + alias])
            command(['docker', 'push', IMAGE + ':' + alias])
        if remote_release['draft']:
            command(['gh', 'release', 'edit', a.tag, '--repo', REPO, '--draft=false', '--latest'])
        state.update(status='COMPLETE', release_url=f'https://github.com/{REPO}/releases/tag/{a.tag}',
                     build_info=info)
        stage('COMPLETE')
    except BaseException as e:
        # Only sanitized exception messages enter durable records.
        state.update(status='FAILED', error=str(e))
        stage('STOPPED_WITHOUT_ROLLBACK')
        raise


if __name__ == '__main__':
    main()
