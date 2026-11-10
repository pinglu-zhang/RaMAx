import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


contract = module('release_contract', ROOT / 'ci/check_release.py')
smoke = module('release_smoke', ROOT / 'ci/release_smoke.py')
glibc = module('check_glibc', ROOT / 'ci/check_glibc.py')


class ReleaseChecks(unittest.TestCase):
    def test_glibc_requires_not_exported_versions(self):
        text = ("Version definition section '.gnu.version_d'\n"
                "  Name: GLIBC_2.39\n"
                "Version needs section '.gnu.version_r'\n"
                "  Name: GLIBC_2.17\n  Name: GLIBC_2.2.5\n  Name: GLIBC_PRIVATE\n")
        self.assertEqual(glibc.required_versions(text), ['GLIBC_2.2.5', 'GLIBC_2.17'])
        self.assertEqual(glibc.required_versions("No version information found."), [])

    def test_current_metadata_and_wrong_tag(self):
        self.assertEqual(contract.read_contract(ROOT, 'v1.1.0'), '1.1.0')
        with self.assertRaises(ValueError):
            contract.read_contract(ROOT, 'v1.0.9')

    def test_metadata_mismatch_and_legacy_helper(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p / 'recipe').mkdir()
            (p / 'CMakeLists.txt').write_text('project(RaMAx VERSION 1.1.0 LANGUAGES CXX)')
            (p / 'Dockerfile').write_text('ARG RAMAX_VERSION=1.1.0\n')
            (p / 'recipe/meta.yaml').write_text('{% set version = "1.0.9" %}\n')
            with self.assertRaises(ValueError):
                contract.read_contract(p)
            (p / 'recipe/meta.yaml').write_text('{% set version = "1.1.0" %}\nhalAppendCactusSubtree\n')
            with self.assertRaises(ValueError):
                contract.read_contract(p)

    def test_full_reverse_walk_and_damaged_paths(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'graph.gfa'
            base = 'H\tVN:Z:1.1\nS\t1\tAC\nS\t2\tGT\nL\t1\t+\t2\t+\t0M\n'
            walks = 'W\tR\t0\tchr\t0\t4\t>1>2\nW\tQ\t0\tchr\t0\t4\t<2<1\n'
            p.write_text(base + walks)
            result = smoke.check_gfa(p, {('R', 'chr'): 'ACGT', ('Q', 'chr'): 'ACGT'})
            self.assertEqual(result['walks'], 2)
            for damaged in [base + walks.replace('<2<1', '<9<1'),
                            base.replace('L\t1\t+\t2\t+\t0M\n', '') + walks,
                            base + walks.replace('\t0\t4\t>1>2', '\t1\t5\t>1>2'),
                            base + walks + walks.splitlines(True)[0]]:
                p.write_text(damaged)
                with self.assertRaises(RuntimeError):
                    smoke.check_gfa(p, {('R', 'chr'): 'ACGT', ('Q', 'chr'): 'ACGT'})

    def test_paf_cigar_and_identity(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'graph.paf'
            row = 'Q.chr\t4\t0\t4\t+\tR.chr\t4\t0\t4\t4\t4\t60\tcg:Z:4=\n'
            p.write_text(row)
            dna = {('Q', 'chr'): 'ACGT', ('R', 'chr'): 'ACGT'}
            self.assertEqual(smoke.check_paf(p, dna)['records'], 1)
            for damaged in [row.replace('4=', '3='), row.replace('Q.chr', 'UNKNOWN.chr'), '']:
                p.write_text(damaged)
                with self.assertRaises(RuntimeError):
                    smoke.check_paf(p, dna)

    def test_maf_row_sequence_and_nonempty_requirement(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'graph.maf'
            dna = 'ACGT' * 4000
            rows = 'a\n' + ''.join(f's {s}.chr 0 16000 + 16000 {dna}\n' for s in ('R', 'A', 'B'))
            expected = {(s, 'chr'): dna for s in ('R', 'A', 'B')}
            p.write_text(rows)
            self.assertEqual(smoke.check_maf(p, expected)['leaf_rows'], 3)
            for damaged in [rows.replace(dna, 'T' + dna[1:], 1), '##maf version=1\n']:
                p.write_text(damaged)
                with self.assertRaises(RuntimeError):
                    smoke.check_maf(p, expected)


if __name__ == '__main__':
    unittest.main()
