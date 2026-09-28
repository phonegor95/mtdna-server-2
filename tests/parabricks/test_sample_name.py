import importlib.util
from pathlib import Path
import unittest
spec = importlib.util.spec_from_file_location('caller', Path(__file__).parents[2] / 'bin/parabricks_mutect.py')
caller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(caller)
class SampleTests(unittest.TestCase):
    def test_same_sample_multiple_groups(self):
        self.assertEqual(caller.sample_name('@RG\tID:a\tSM:sample\n@RG\tID:b\tSM:sample\n'), 'sample')
    def test_reject_ambiguous_or_absent_samples(self):
        for header in ('@SQ\tSN:chrM', '@RG\tID:a', '@RG\tID:a\tSM:', '@RG\tID:a\tSM:a\n@RG\tID:b\tSM:b'):
            with self.subTest(header=header), self.assertRaises(ValueError):
                caller.sample_name(header)
if __name__ == '__main__': unittest.main()
