from .create_data import make_random_genome, make_contigs, \
     write_genome_and_contigs
from .utils import TempStuff, runquast, load_report, assert_match

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
    
