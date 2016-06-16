import tempfile
import os.path
import sys
from cStringIO import StringIO

p = os.path.dirname(__file__)
p = os.path.join(p, '..')

sys.path.insert(0, p)
import quast

class TempStuff(object):
    def __init__(self, suff=''):
        self.dirname = tempfile.mkdtemp(suffix=suff)

    def filename(self, name):
        path = os.path.join(self.dirname, name)
        return path

def runquast(args):
    return quast.main(args)
        

def load_report(output_dir):
    reportfile = os.path.join(output_dir, 'report.tsv')
    d = {}
    
    firstline = True
    for line in open(reportfile):
        if firstline:
            firstline = False
            continue
        line = line.strip().split('\t', 1)
        try:
            d[line[0]] = float(line[1])
        except ValueError:
            d[line[0]] = line[1]
    return d

def assert_match(report, field, value):
    actual = report.get(field)
    assert value == actual, (field, value, actual)
