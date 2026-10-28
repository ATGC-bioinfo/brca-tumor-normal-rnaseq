#!/usr/bin/env python3
"""Convert a SnpEff (+ SnpSift ClinVar) annotated Mutect2 VCF into a MAF file for maftools.

One MAF row per variant (first ALT allele). The SnpEff annotation used is the first
transcript-level ANN entry for that allele (SnpEff lists the most severe effect first).
Coordinates follow the MAF convention: indels drop the shared anchor base, '-' marks an empty allele.

Usage: snpeff_to_maf.py in.vcf.gz --tumor T1 --normal N1 --barcode P1 > out.maf
"""
import argparse
import gzip
import sys

# SnpEff / Sequence Ontology term -> MAF Variant_Classification (first matching term wins)
SO_TO_MAF = {
    'stop_gained': 'Nonsense_Mutation',
    'stop_lost': 'Nonstop_Mutation',
    'start_lost': 'Translation_Start_Site',
    'initiator_codon_variant': 'Translation_Start_Site',
    'splice_acceptor_variant': 'Splice_Site',
    'splice_donor_variant': 'Splice_Site',
    'frameshift_variant': 'Frame_Shift_',          # + Del / Ins
    'inframe_deletion': 'In_Frame_Del',
    'disruptive_inframe_deletion': 'In_Frame_Del',
    'conservative_inframe_deletion': 'In_Frame_Del',
    'inframe_insertion': 'In_Frame_Ins',
    'disruptive_inframe_insertion': 'In_Frame_Ins',
    'conservative_inframe_insertion': 'In_Frame_Ins',
    'missense_variant': 'Missense_Mutation',
    'splice_region_variant': 'Splice_Region',
    'splice_donor_region_variant': 'Splice_Region',
    'splice_donor_5th_base_variant': 'Splice_Region',
    'splice_polypyrimidine_tract_variant': 'Splice_Region',
    'synonymous_variant': 'Silent',
    'stop_retained_variant': 'Silent',
    'start_retained_variant': 'Silent',
    '5_prime_UTR_variant': "5'UTR",
    '5_prime_UTR_premature_start_codon_gain_variant': "5'UTR",
    '3_prime_UTR_variant': "3'UTR",
    'non_coding_transcript_exon_variant': 'RNA',
    'non_coding_transcript_variant': 'RNA',
    'intron_variant': 'Intron',
    'upstream_gene_variant': "5'Flank",
    'downstream_gene_variant': "3'Flank",
    'intergenic_region': 'IGR',
}

AA3 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E',
       'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
       'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*',
       'Sec': 'U', 'Xaa': 'X'}

COLS = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_Position', 'End_Position',
        'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele',
        'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode',
        'Matched_Norm_Sample_Barcode', 'HGVSc', 'HGVSp_Short', 'Transcript_ID', 'IMPACT',
        'SnpEff_Effect', 't_depth', 't_ref_count', 't_alt_count', 't_vaf',
        'n_depth', 'n_ref_count', 'n_alt_count', 'ClinVar_CLNSIG', 'ClinVar_CLNDN',
        'ClinVar_CLNREVSTAT', 'ClinVar_ONC', 'ClinVar_ONCDN', 'ClinVar_SCI', 'ClinVar_SCIDN']


def hgvsp_short(p):
    """p.His1047Arg -> p.H1047R"""
    for three, one in AA3.items():
        p = p.replace(three, one)
    return p


def maf_alleles(pos, ref, alt):
    """VCF (pos, ref, alt) -> (start, end, ref, alt, type) in MAF convention."""
    if len(ref) == len(alt):
        vtype = {1: 'SNP', 2: 'DNP', 3: 'TNP'}.get(len(ref), 'ONP')
        return pos, pos + len(ref) - 1, ref, alt, vtype
    if ref[0] == alt[0]:  # drop the shared anchor base
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if len(ref) > len(alt):
        return pos, pos + len(ref) - 1, ref, alt or '-', 'DEL'
    if not ref:  # pure insertion: MAF start/end are the bases flanking it
        return pos - 1, pos, '-', alt, 'INS'
    return pos, pos + len(ref) - 1, ref, alt, 'INS'


def classify(effect, vtype):
    for term in effect.split('&'):
        if term in SO_TO_MAF:
            c = SO_TO_MAF[term]
            return c + ('Del' if vtype == 'DEL' else 'Ins') if c == 'Frame_Shift_' else c
    return 'Targeted_Region'


def pick_ann(info, alt):
    anns = [a.split('|') for a in info.get('ANN', '').split(',') if a]
    anns = [a for a in anns if len(a) > 10 and a[0] == alt]
    tx = [a for a in anns if a[5] == 'transcript']
    return (tx or anns or [None])[0]


def convert(lines, tumor, normal, barcode):
    """Yield MAF rows (dicts) from VCF lines."""
    ti = ni = None
    for line in lines:
        if line.startswith('##'):
            continue
        f = line.rstrip('\n').split('\t')
        if line.startswith('#'):
            ti, ni = f.index(tumor), f.index(normal)
            continue
        chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4].split(',')[0]
        info = dict(x.split('=', 1) if '=' in x else (x, '') for x in f[7].split(';'))
        fmt = f[8].split(':')
        t = dict(zip(fmt, f[ti].split(':')))
        n = dict(zip(fmt, f[ni].split(':')))
        t_ad, n_ad = t.get('AD', '.,.').split(','), n.get('AD', '.,.').split(',')
        start, end, mref, malt, vtype = maf_alleles(pos, ref, alt)
        ann = pick_ann(info, alt)
        effect = ann[1] if ann else 'intergenic_region'
        yield {
            'Hugo_Symbol': (ann[3] if ann and ann[3] else 'Unknown'),
            'Entrez_Gene_Id': '0', 'Chromosome': chrom,
            'Start_Position': start, 'End_Position': end, 'Strand': '+',
            'Variant_Classification': classify(effect, vtype), 'Variant_Type': vtype,
            'Reference_Allele': mref, 'Tumor_Seq_Allele1': mref, 'Tumor_Seq_Allele2': malt,
            'Tumor_Sample_Barcode': barcode, 'Matched_Norm_Sample_Barcode': normal,
            'HGVSc': ann[9] if ann else '', 'HGVSp_Short': hgvsp_short(ann[10]) if ann else '',
            'Transcript_ID': ann[6] if ann else '', 'IMPACT': ann[2] if ann else 'MODIFIER',
            'SnpEff_Effect': effect,
            't_depth': t.get('DP', ''), 't_ref_count': t_ad[0], 't_alt_count': t_ad[1],
            't_vaf': t.get('AF', '').split(',')[0],
            'n_depth': n.get('DP', ''), 'n_ref_count': n_ad[0], 'n_alt_count': n_ad[1],
            'ClinVar_CLNSIG': info.get('CLNSIG', ''), 'ClinVar_CLNDN': info.get('CLNDN', ''),
            'ClinVar_CLNREVSTAT': info.get('CLNREVSTAT', ''),
            'ClinVar_ONC': info.get('ONC', ''), 'ClinVar_ONCDN': info.get('ONCDN', ''),
            'ClinVar_SCI': info.get('SCI', ''), 'ClinVar_SCIDN': info.get('SCIDN', ''),
        }


def demo():
    assert maf_alleles(100, 'A', 'G') == (100, 100, 'A', 'G', 'SNP')
    assert maf_alleles(100, 'AT', 'A') == (101, 101, 'T', '-', 'DEL')
    assert maf_alleles(100, 'A', 'AGG') == (100, 101, '-', 'GG', 'INS')
    assert maf_alleles(100, 'AC', 'GT') == (100, 101, 'AC', 'GT', 'DNP')
    assert classify('frameshift_variant&splice_region_variant', 'DEL') == 'Frame_Shift_Del'
    assert classify('missense_variant&splice_region_variant', 'SNP') == 'Missense_Mutation'
    assert classify('something_new', 'SNP') == 'Targeted_Region'
    assert hgvsp_short('p.His1047Arg') == 'p.H1047R'
    assert hgvsp_short('p.Arg213Ter') == 'p.R213*'
    vcf = ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN1\tT1\n',
           'chr3\t179234297\t.\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;ANN=G|missense_variant|MODERATE|PIK3CA|E|transcript|ENST1|protein_coding|21/21|c.3140A>G|p.His1047Arg|||||\tGT:AD:AF:DP\t0/0:30,0:0.01:30\t0/1:40,10:0.2:50\n']
    r = next(convert(vcf, 'T1', 'N1', 'P1'))
    assert (r['Hugo_Symbol'], r['HGVSp_Short'], r['t_alt_count'], r['ClinVar_CLNSIG']) == ('PIK3CA', 'p.H1047R', '10', 'Pathogenic')
    print('snpeff_to_maf demo OK')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('vcf', nargs='?')
    p.add_argument('--tumor')
    p.add_argument('--normal')
    p.add_argument('--barcode')
    p.add_argument('--demo', action='store_true')
    a = p.parse_args()
    if a.demo:
        return demo()
    out = sys.stdout
    out.write('\t'.join(COLS) + '\n')
    with gzip.open(a.vcf, 'rt') as fh:
        for row in convert(fh, a.tumor, a.normal, a.barcode):
            out.write('\t'.join(str(row[c]) for c in COLS) + '\n')


if __name__ == '__main__':
    main()
