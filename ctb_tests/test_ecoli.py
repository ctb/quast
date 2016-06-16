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

    is_diff = False
    old_results_subdir = 'results-ecoli-5m-dn/'
    for filename in [
        'contigs_reports/alignments_megahit_5m_dn.tsv',
        'contigs_reports/contigs_report_megahit_5m_dn.mis_contigs.info',
        'contigs_reports/megahit_5m_dn.mis_contigs.fa',
        'contigs_reports/misassemblies_report.tsv',
        'contigs_reports/misassemblies_report.txt',
        'contigs_reports/nucmer_output/megahit_5m_dn.all_snps',
        'contigs_reports/nucmer_output/megahit_5m_dn.coords',
        'contigs_reports/nucmer_output/megahit_5m_dn.coords.filtered',
        'contigs_reports/nucmer_output/megahit_5m_dn.unaligned',
        'contigs_reports/nucmer_output/megahit_5m_dn.used_snps',
        'contigs_reports/transposed_report_misassemblies.tsv',
        'contigs_reports/transposed_report_misassemblies.txt',
        'contigs_reports/unaligned_report.tsv',
        'contigs_reports/unaligned_report.txt',
        'genome_stats/genome_info.txt',
        'genome_stats/megahit_5m_dn_gaps.txt',
        'report.tsv',
        'report.txt',
        'transposed_report.tsv',
        'transposed_report.txt',]:
        new_file = tmp.filename(os.path.join('output', filename))
        old_file = get_test_data(os.path.join(old_results_subdir, filename))

        
        if open(new_file).read() != open(old_file).read():
            print '%s differs from %s' % (new_file, old_file)
            is_diff = True
    if is_diff:
        assert 0, 'differences exist'


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
    
