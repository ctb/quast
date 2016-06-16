from .utils import TempStuff, runquast, load_report, assert_match
import pprint
import screed
import os.path


#ecoliMG1655.fa.gz
#megahit_1m.fa.gz
#megahit_500k.fa.gz
#megahit_5m.fa.gz
#megahit_5m_dn.fa.gz


def get_test_data(filename):
    thisdir = os.path.dirname(os.path.normpath(__file__))
    datadir = os.path.join(thisdir, 'test-data')
    fullpath = os.path.join(datadir, filename)
    return fullpath


def test_1():
    tmp = TempStuff('real')
    reference = get_test_data('ecoliMG1655.fa.gz')
    contigs = get_test_data('megahit_5m_dn.fa.gz')

    runquast(['-R', reference,
              '-o', tmp.filename('output'),
              contigs])


def test_2():
    tmp = TempStuff('real')
    reference = get_test_data('ecoliMG1655.fa.gz')
    contigs = get_test_data('megahit_1m.fa.gz')

    runquast(['-R', reference,
              '-o', tmp.filename('output'),
              contigs])


def test_3():
    tmp = TempStuff('real')
    reference = get_test_data('ecoliMG1655.fa.gz')
    contigs = get_test_data('megahit_500k.fa.gz')

    runquast(['-R', reference,
              '-o', tmp.filename('output'),
              contigs])


def test_4():
    tmp = TempStuff('real')
    reference = get_test_data('ecoliMG1655.fa.gz')
    contigs = get_test_data('megahit_5m.fa.gz')

    runquast(['-R', reference,
              '-o', tmp.filename('output'),
              contigs])
    
