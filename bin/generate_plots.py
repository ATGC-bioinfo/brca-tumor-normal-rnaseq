#!/usr/bin/env python3
"""Interpretation plots for the somatic calls (one PNG + PDF each).

1. filter_tumor_vs_normal_vaf  every Mutect2 call, tumor VAF vs normal VAF, coloured by final status:
                               shows the filters separating somatic calls from germline / noise.
2. mrna_strand_spectrum        substitutions written on the gene's mRNA strand. RNA editing (ADAR)
                               shows up as mRNA A>G far above T>C; DNA mutations have no such imbalance.
3. vaf_distribution            tumor VAF density of final calls per patient: a peak near the tumor's
                               clonal VAF = truncal mutations, a low-VAF tail = subclones.

Usage: generate_plots.py --filtered P*.somatic.filtered.vcf.gz --pass P*.somatic.pass.vcf.gz --gtf genes.gtf
"""
import argparse
import bisect
import gzip
import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

COMP = str.maketrans('ACGT', 'TGCA')
PAIRS = [('A>G', 'T>C'), ('C>T', 'G>A'), ('C>A', 'G>T'), ('A>C', 'T>G'), ('A>T', 'T>A'), ('C>G', 'G>C')]


def read_calls(path):
    """Yield (chrom, pos, ref, alt, filter, tumor_af, normal_af) from a Mutect2 VCF."""
    tumor = None
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##tumor_sample='):
                tumor = line.strip().split('=', 1)[1]
            elif line.startswith('#CHROM'):
                cols = line.rstrip('\n').split('\t')
                ti = cols.index(tumor)
                ni = 9 if ti == 10 else 10
            elif not line.startswith('#'):
                f = line.rstrip('\n').split('\t')
                fmt = f[8].split(':')
                t, n = dict(zip(fmt, f[ti].split(':'))), dict(zip(fmt, f[ni].split(':')))
                try:
                    taf, naf = float(t['AF'].split(',')[0]), float(n['AF'].split(',')[0])
                except (KeyError, ValueError):
                    continue
                yield f[0], int(f[1]), f[3], f[4].split(',')[0], f[6], taf, naf


def gene_strands(gtf):
    """chrom -> (sorted starts, [(start, end, strand)]) for protein-coding genes."""
    genes = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.split('\t', 8)
            if f[2] == 'gene' and 'gene_type "protein_coding"' in f[8]:
                genes[f[0]].append((int(f[3]), int(f[4]), f[6]))
    out = {}
    for c, g in genes.items():
        g.sort()
        out[c] = ([s for s, _, _ in g], g)
    return out


def strand_at(index, chrom, pos, max_gene=2_500_000):
    """Strand of the protein-coding gene(s) covering pos, or None if none / both strands."""
    if chrom not in index:
        return None
    starts, genes = index[chrom]
    i = bisect.bisect_right(starts, pos)
    s = {st for a, b, st in genes[bisect.bisect_left(starts, pos - max_gene):i] if a <= pos <= b}
    return s.pop() if len(s) == 1 else None


def mrna_change(ref, alt, strand):
    """Substitution as read on the mRNA: unchanged for + genes, complemented for - genes."""
    return f'{ref}>{alt}' if strand == '+' else f'{ref.translate(COMP)}>{alt.translate(COMP)}'


def save(fig, outdir, name):
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(outdir, f'{name}.{ext}'), dpi=200)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--filtered', nargs='+', required=True)
    p.add_argument('--pass', dest='passed', nargs='+', required=True)
    p.add_argument('--gtf', required=True)
    p.add_argument('--outdir', default='.')
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    pid = lambda path: os.path.basename(path).split('.')[0]  # noqa: E731
    filtered = sorted(a.filtered, key=pid)
    passed = {pid(v): list(read_calls(v)) for v in a.passed}
    patients = sorted(passed)

    # 1. Tumor vs normal VAF, all calls
    fig, axes = plt.subplots(1, len(filtered), figsize=(4.2 * len(filtered), 4.2), squeeze=False, sharey=True)
    for ax, vcf in zip(axes[0], filtered):
        final = {c[:4] for c in passed.get(pid(vcf), [])}
        groups = {'failed GATK filters': [], 'GATK PASS, failed project cutoffs': [], 'final somatic call': []}
        for c in read_calls(vcf):
            key = 'final somatic call' if c[:4] in final else \
                  'GATK PASS, failed project cutoffs' if c[4] == 'PASS' else 'failed GATK filters'
            groups[key].append((c[6], c[5]))
        for (label, pts), col, size in zip(groups.items(), ('#c8c8c8', '#f0a33a', '#c0392b'), (3, 6, 10)):
            if pts:
                ax.scatter(*zip(*pts), s=size, c=col, alpha=0.6, lw=0, label=f'{label} ({len(pts)})')
        ax.set(title=pid(vcf), xlabel='Normal VAF', xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
        ax.legend(fontsize=7, loc='lower right', markerscale=2)
    axes[0][0].set_ylabel('Tumor VAF')
    fig.suptitle('Every Mutect2 call: final somatic calls (red) are tumor-only and well supported; orange = removed by project cutoffs')
    save(fig, a.outdir, 'filter_tumor_vs_normal_vaf')

    # 2. mRNA-strand spectrum of final SNVs
    index = gene_strands(a.gtf)
    counts = {}
    for pt in patients:
        c = Counter()
        for chrom, pos, ref, alt, *_ in passed[pt]:
            st = strand_at(index, chrom, pos) if len(ref) == len(alt) == 1 else None
            if st:
                c[mrna_change(ref, alt, st)] += 1
        counts[pt] = c
    fig, ax = plt.subplots(figsize=(10, 4.5))
    width = 0.8 / max(len(patients), 1)
    for j, pt in enumerate(patients):
        for i, (x, y) in enumerate(PAIRS):
            ax.bar(i * 2 + j * width, counts[pt][x], width, color=plt.cm.tab10(j), label=pt if i == 0 else None)
            ax.bar(i * 2 + 1 + j * width, counts[pt][y], width, color=plt.cm.tab10(j), alpha=0.45)
    ax.set_xticks([k + 0.4 - width / 2 for k in range(12)], [s for pr in PAIRS for s in pr])
    ax.set(ylabel='SNVs in protein-coding genes', xlabel='Substitution on the mRNA strand (each pair = complementary changes)')
    ax.set_title('RNA editing check: ADAR editing would make mRNA A>G (first bar) much taller than T>C')
    ax.legend()
    save(fig, a.outdir, 'mrna_strand_spectrum')

    # 3. Tumor VAF distribution of final calls
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for j, pt in enumerate(patients):
        vafs = [c[5] for c in passed[pt]]
        ax.hist(vafs, bins=[i / 25 for i in range(26)], histtype='step', lw=2, color=plt.cm.tab10(j),
                label=f'{pt} (n={len(vafs)}, median {sorted(vafs)[len(vafs) // 2]:.2f})' if vafs else pt)
    ax.set(xlabel='Tumor VAF', ylabel='Somatic calls', xlim=(0, 1),
           title='Clonal mutations cluster at the tumor\'s purity/2; a low-VAF tail = subclones')
    ax.legend()
    save(fig, a.outdir, 'vaf_distribution')


def demo():
    assert mrna_change('A', 'G', '+') == 'A>G' and mrna_change('T', 'C', '-') == 'A>G'
    idx = {'chr1': ([100, 500], [(100, 400, '+'), (500, 900, '-')])}
    assert strand_at(idx, 'chr1', 150) == '+' and strand_at(idx, 'chr1', 600) == '-'
    assert strand_at(idx, 'chr1', 450) is None and strand_at(idx, 'chr2', 1) is None
    print('generate_plots demo OK')


if __name__ == '__main__':
    import sys
    demo() if '--demo' in sys.argv else main()
