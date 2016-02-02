from .create_data import make_random_genome, make_contigs, \
     write_genome_and_contigs
from .utils import TempStuff, runquast, load_report, assert_match
import pprint

def test_simple():
    L2 = 10000
    genome = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome, contig_ivals)

    tmp = TempStuff('test1')
    write_genome_and_contigs(tmp.dirname, [genome], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    assert_match(report, 'Largest alignment', 10000)

    #outfilename = 'test1_out/contigs_reports/contigs_report_contigs.stdout'
    #out = open(outfilename).read()
    #start = out.find('Uncovered Regions')
    #print out[start:start+55]
    
    
def test_simple_2():
    L2 = 10000

    genome1 = make_random_genome(1e5)
    genome2 = make_random_genome(1e5)
    genome3 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)
    contigs.extend(make_contigs(genome2, contig_ivals))

    # make a contig that is erroneous --
    last_contig = contigs[-1]
    contigs = contigs[:-1]
    last_contig = last_contig[L2/2:] + last_contig[:L2/2]
    contigs.append(last_contig)

    # make a contig that aligns to nothing at all
    noalign = make_random_genome(L2)
    contigs.append(noalign)

    # make a contig part of which aligns to nothing.
    noalign = make_random_genome(L2/2)
    contigs[0] = noalign + contigs[0][:L2/2]

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1, genome2, genome3], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, 'Largest alignment', 10000)
    assert_match(report, '# unaligned contigs', '1 + 1 part')

    #outfilename = 'test2_out/contigs_reports/contigs_report_contigs.stdout'
    #out = open(outfilename).read()
    #start = out.find('Uncovered Regions')
    #print out[start:start+55]
