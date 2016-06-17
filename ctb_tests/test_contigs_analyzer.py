from libs.contigs_analyzer import plantakolya_analyze_coverage
import StringIO
import recordtype

Alignment = recordtype.recordtype('alignmentClass',
                                  ['s1', 's2', 'e1', 'e2', 'contig'])


def test_empty():
    regions = {}
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None
    
    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    assert nothing_aligned
    # non-error execution is good enough for me...


def test_empty_ref_aligns():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': None }
    ref_aligns = {}

    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    assert 'ERROR: Reference someref does not have any alignments!' in logtext


def test_2():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': [(0, 100)] }
    ref_aligns = {}
    ref_aligns['someref'] = [ Alignment(s1=0, s2=0, e1=100, e2=100,
                                        contig='contig1') ]


    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    print logtext
    assert 'Reference someref: 1 total alignments. 1 total regions.' in logtext


def test_3():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': [(75, 100)] }
    ref_aligns = {}
    ref_aligns['someref'] = [ Alignment(s1=0, s2=0, e1=50, e2=50,
                                        contig='contig1') ]


    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    print logtext
    assert 'Reference someref: 1 total alignments. 1 total regions.' in logtext
    assert 'This align occurs before our region of interest, skip' in logtext
    assert 'There are no more aligns. Skipping this region.' in logtext


def test_4():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': [(25, 100)] }
    ref_aligns = {}
    ref_aligns['someref'] = [ Alignment(s1=0, s2=0, e1=50, e2=50,
                                        contig='contig1') ]


    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    print logtext
    assert 'Reference someref: 1 total alignments. 1 total regions.' in logtext
    assert 'START within alignment' in logtext
    assert 'Moving reference start from 0 to 25' in logtext
    assert 'Moving contig start from 0 to 25' in logtext
    # assert 0 is this right??


def test_5():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': [(0, 25)] }
    ref_aligns = {}
    ref_aligns['someref'] = [ Alignment(s1=50, s2=50, e1=100, e2=100,
                                        contig='contig1') ]


    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    print logtext
    assert 'Reference someref: 1 total alignments. 1 total regions.' in logtext
    assert 'There are no aligns within this region.' in logtext
    # assert 0 is this right??


def test_6():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': [(0, 25)] }
    ref_aligns = {}
    ref_aligns['someref'] = [ Alignment(s1=10, s2=10, e1=100, e2=100,
                                        contig='contig1') ]


    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    print logtext
    assert 'Reference someref: 1 total alignments. 1 total regions.' in logtext
    assert 'START in gap: 0 to 10 (10 bp)' in logtext
    assert 'END in current alignment.  Modifying 100 to 25.' in logtext
    # assert 0 is this right??


def test_7():
    planta_out_f = StringIO.StringIO()
    total_indels_info = None
    ref_aligns = {}
    ref_features = {}
    snps = []
    contig_ns = []
    used_snps_file = None

    regions = { 'someref': [(25, 100)] }
    ref_aligns = {}
    ref_aligns['someref'] = [ Alignment(s1=0, s2=50, e1=50, e2=0,
                                        contig='contig1') ]


    total_indels_info, region_covered, region_ambig, gaps, neg_gaps, \
           redundant, nothing_aligned, total_redundant = \
       plantakolya_analyze_coverage(regions, planta_out_f, total_indels_info,
                                    ref_aligns, ref_features, snps, contig_ns,
                                    used_snps_file)

    logtext = planta_out_f.getvalue()
    print logtext
    assert 'Reference someref: 1 total alignments. 1 total regions.' in logtext
    assert 'START within alignment' in logtext
    assert 'Moving reference start from 0 to 25' in logtext
    assert 'Moving contig start from 50 to 25' in logtext
