import importlib.util
import os
from tempfile import NamedTemporaryFile


# Load visualization module from file path because of hyphen in directory name
module_path = os.path.join(os.path.dirname(__file__), '..', 'AviTag-Seq_0411', 'visualization.py')
spec = importlib.util.spec_from_file_location('visualization', module_path)
visualization = importlib.util.module_from_spec(spec)
spec.loader.exec_module(visualization)


def test_parse_sites_file():
    """Ensure parseSitesFile returns sorted results and reference sequence."""
    header = '\t'.join(['h'] * 33) + '\n'
    # two rows with the same sequence to test aggregation
    row1 = ['chr1:100', 'x','x','x','x','x','x','x','x','x','x','5'] + ['x']*9 + ['AAAA'] + ['x']*10 + ['AAAANAAA']
    row2 = ['chr1:100', 'x','x','x','x','x','x','x','x','x','x','3'] + ['x']*9 + ['AAAA'] + ['x']*10 + ['AAAANAAA']
    line1 = '\t'.join(row1) + '\n'
    line2 = '\t'.join(row2) + '\n'
    with NamedTemporaryFile('w+', delete=False) as tmp:
        tmp.write(header)
        tmp.write(line1)
        tmp.write(line2)
        tmp_path = tmp.name
    try:
        offtargets, ref_seq, chrs = visualization.parseSitesFile(tmp_path)
        assert ref_seq.strip() == 'AAAANAAA'
        assert offtargets[0]['reads'] == 8
        assert chrs[0].startswith('1:')
    finally:
        os.remove(tmp_path)
