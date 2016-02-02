import tempfile
import os.path
import sys
from cStringIO import StringIO
import nose

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
