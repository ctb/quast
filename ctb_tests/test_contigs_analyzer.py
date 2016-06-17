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
    assert 0
