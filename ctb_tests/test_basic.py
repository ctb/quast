from .create_data import make_random_genome, make_contigs, \
     write_genome_and_contigs
from .utils import TempStuff, runquast, load_report, assert_match
import pprint
import screed

def test_basic():
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
    assert_match(report, '# unaligned contigs', '0 + 0 part')
    assert_match(report, '# contigs (>= 10000 bp)', 2)
    assert_match(report, '# contigs (>= 50000 bp)', 0)
    assert_match(report, '# contigs', 2)
    assert_match(report, 'Duplication ratio', 1.0)
    assert_match(report, 'Genome fraction (%)', 20.0)
    assert_match(report, 'Reference GC (%)', 50.0)
    assert_match(report, 'Reference length', 100000.0)
    assert_match(report, 'Largest contig', 10000.0)

    #outfilename = 'test1_out/contigs_reports/contigs_report_contigs.stdout'
    #out = open(outfilename).read()
    #start = out.find('Uncovered Regions')
    #print out[start:start+55]


def test_basic_rc():
    L2 = 10000
    genome = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome, contig_ivals)

    contigs[0] = screed.rc(contigs[0])

    tmp = TempStuff('test1')
    write_genome_and_contigs(tmp.dirname, [genome], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    assert_match(report, 'Largest alignment', 10000)
    assert_match(report, '# unaligned contigs', '0 + 0 part')
    assert_match(report, '# contigs (>= 10000 bp)', 2)
    assert_match(report, '# contigs (>= 50000 bp)', 0)
    assert_match(report, '# contigs', 2)
    assert_match(report, 'Duplication ratio', 1.0)
    assert_match(report, 'Genome fraction (%)', 20.0)
    assert_match(report, 'Reference GC (%)', 50.0)
    assert_match(report, 'Reference length', 100000.0)
    assert_match(report, 'Largest contig', 10000.0)

    #outfilename = 'test1_out/contigs_reports/contigs_report_contigs.stdout'
    #out = open(outfilename).read()
    #start = out.find('Uncovered Regions')
    #print out[start:start+55]
    
    
def test_largest_contig():
    L2 = 10000
    genome = make_random_genome(1e5)

    contig_ivals = [(10000, L2 + 5), (50000, L2)]
    contigs = make_contigs(genome, contig_ivals)

    tmp = TempStuff('test1')
    write_genome_and_contigs(tmp.dirname, [genome], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    assert_match(report, 'Largest alignment', 10005)
    assert_match(report, '# contigs', 2)
    assert_match(report, 'Largest contig', 10005.0)

    
def test_multiple_genomes():
    L2 = 10000

    genome1 = make_random_genome(1e5)
    genome2 = make_random_genome(1e5)
    genome3 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)
    contigs.extend(make_contigs(genome2, contig_ivals))

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1, genome2, genome3], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, 'Largest alignment', 10000)
    assert_match(report, '# unaligned contigs', '0 + 0 part')
    assert_match(report, '# contigs (>= 10000 bp)', 4)
    assert_match(report, '# contigs (>= 50000 bp)', 0)
    assert_match(report, '# contigs', 4)
    assert_match(report, 'Reference GC (%)', 50.0)
    assert_match(report, 'Reference length', 300000.0)


def test_noalign():
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)

    # make a contig that aligns to nothing at all
    noalign = make_random_genome(L2)
    contigs.append(noalign)

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, 'Largest alignment', 10000)
    assert_match(report, '# unaligned contigs', '1 + 0 part')
    assert_match(report, '# contigs (>= 10000 bp)', 3)
    assert_match(report, '# contigs (>= 50000 bp)', 0)
    assert_match(report, '# contigs', 3)


def test_misassembly():
    import random
    random.seed(1)                      # otherwise, "Largest alignment" varies
    
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)

    # make a misassembled contig
    last_contig = contigs[-1]
    contigs = contigs[:-1]
    last_contig = last_contig[L2/2:] + last_contig[:L2/2]
    contigs.append(last_contig)

    # make a contig part of which aligns to nothing.
    noalign = make_random_genome(L2/2)
    contigs[0] = noalign + contigs[0][:L2/2]

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, 'Largest alignment', 5002)
    assert_match(report, '# misassemblies', 1)
    assert_match(report, '# misassembled contigs', 1.0)
    assert_match(report, 'Misassembled contigs length', 10000.0)
    assert_match(report, '# unaligned contigs', '0 + 1 part')
    assert_match(report, '# contigs (>= 10000 bp)', 2)
    assert_match(report, '# contigs (>= 50000 bp)', 0)
    assert_match(report, '# contigs', 2)
    assert_match(report, 'Unaligned length', 5000)
    assert_match(report, 'Total length', 20000.0)
    assert_match(report, 'Total length (>= 0 bp)', 20000.0)
    assert_match(report, 'Total length (>= 1000 bp)', 20000.0)
    assert_match(report, 'Total length (>= 10000 bp)', 20000.0)
    assert_match(report, 'Total length (>= 25000 bp)', 0.0)
    assert_match(report, 'Total length (>= 5000 bp)', 20000.0)
    assert_match(report, 'Total length (>= 50000 bp)', 0.0)
    assert_match(report, '# contigs (>= 0 bp)', 2.0)
    assert_match(report, '# contigs (>= 1000 bp)', 2.0)
    assert_match(report, '# contigs (>= 25000 bp)', 0.0)
    assert_match(report, '# contigs (>= 5000 bp)', 2.0)
    assert_match(report, 'Genome fraction (%)', 15.0)


def test_misassembly_plus_unaligned():
    import random
    random.seed(1)                      # otherwise, "Largest alignment" varies

    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)

    # make a misassembled contig
    last_contig = contigs[-1]
    contigs = contigs[:-1]
    last_contig = last_contig[L2/2:] + last_contig[:L2/2]
    contigs.append(last_contig)

    # make a contig part of which aligns to nothing.
    noalign = make_random_genome(L2/2)
    contigs[0] = noalign + contigs[0][:L2/2]

    # make another contig, none of which aligns
    contigs.append(make_random_genome(L2))

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, 'Largest alignment', 5002)
    assert_match(report, '# misassemblies', 1)
    assert_match(report, '# misassembled contigs', 1.0)
    assert_match(report, 'Misassembled contigs length', 10000.0)
    assert_match(report, '# unaligned contigs', '1 + 1 part')
    assert_match(report, '# contigs (>= 10000 bp)', 3)
    assert_match(report, '# contigs (>= 50000 bp)', 0)
    assert_match(report, '# contigs', 3)
    assert_match(report, 'Unaligned length', 15000)
    assert_match(report, 'Total length', 30000.0)
    assert_match(report, 'Total length (>= 0 bp)', 30000.0)
    assert_match(report, 'Total length (>= 1000 bp)', 30000.0)
    assert_match(report, 'Total length (>= 10000 bp)', 30000.0)
    assert_match(report, 'Total length (>= 25000 bp)', 0.0)
    assert_match(report, 'Total length (>= 5000 bp)', 30000.0)
    assert_match(report, 'Total length (>= 50000 bp)', 0.0)
    assert_match(report, '# contigs (>= 0 bp)', 3.0)
    assert_match(report, '# contigs (>= 1000 bp)', 3.0)
    assert_match(report, '# contigs (>= 25000 bp)', 0.0)
    assert_match(report, '# contigs (>= 5000 bp)', 3.0)
    assert_match(report, 'Genome fraction (%)', 15.0)


def test_ns_per_kb():
    import random
    random.seed(2)                      # otherwise, "Largest alignment" varies
    
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)

    contig = list(contigs[0])
    for i in range(1, len(contig), 10):
        contig[i] = 'N'
    contig = "".join(contig)
    contigs[0] = contig

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    num_p100kb = contig.count('N') * float(100000) / float(len(contig)) / \
       len(contigs)
    assert_match(report, '# N\'s per 100 kbp', num_p100kb)


def test_indels_insertions():
    import random
    random.seed(1)                      # otherwise, "Largest alignment" varies
    
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)
    
    contig = list(contigs[0])
    total_indels = 0
    for i in range(L2 - 10, 0, -100):
        contig.insert(i, random.choice(['a', 'c', 'g', 't']))
        total_indels += 1
    contig = "".join(contig)
    contigs[0] = contig

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, '# indels per 100 kbp', 500)


def test_indels_deletions():
    import random
    random.seed(1)                      # otherwise, "Largest alignment" varies
    
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)
    
    contig = list(contigs[0])
    total_indels = 0
    for i in range(L2 - 10, 0, -100):
        contig.pop(i)
        total_indels += 1
    contig = "".join(contig)
    contigs[0] = contig

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, '# indels per 100 kbp', 500)


def test_mismatches():
    import random
    random.seed(1)                      # otherwise, "Largest alignment" varies
    
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)
    
    contig = list(contigs[0])
    total_mismatches = 0
    for i in range(L2 - 10, 0, -100):
        ch = contig[i]
        while ch == contig[i]:
            ch = random.choice(['A', 'C', 'G', 'T'])
        contig[i] = ch
        total_mismatches += 1
    contig = "".join(contig)
    contigs[0] = contig

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    print total_mismatches
    assert_match(report, '# mismatches per 100 kbp', 500)


def test_duplication_ratio():
    import random
    random.seed(1)                      # otherwise, "Largest alignment" varies
    
    L2 = 10000

    genome1 = make_random_genome(1e5)

    contig_ivals = [(10000, L2), (50000, L2)]
    contigs = make_contigs(genome1, contig_ivals)
    contigs.append(contigs[0])

    tmp = TempStuff('test2')
    write_genome_and_contigs(tmp.dirname, [genome1], contigs)

    runquast(['-R', tmp.filename('genome.fa'),
              '-o', tmp.filename('output'),
              tmp.filename('contigs.fa')])

    report = load_report(tmp.filename('output'))
    pprint.pprint(report)
    assert_match(report, 'Duplication ratio', 1.5)


{'# local misassemblies': 0.0,
 'GC (%)': 49.81,
 'L50': 1.0,
 'L75': 2.0,
 'LA50': 2.0,
 'LA75': 3.0,
 'N50': 10000.0,
 'N75': 10000.0,
 'NA50': 5000.0,
 'NA75': 4999.0,
 'NGA50': '-',
}
