#!/usr/bin/env python3
"""Keep only Mutect2 calls that are FILTER=PASS and meet our depth / VAF cutoffs.

Reads a FilterMutectCalls VCF (.vcf.gz), writes plain VCF to stdout.
Multi-allelic sites: the strongest alt allele is used for each check.
A missing value ('.') fails the variant.
Optional --editing_sites (REDIportal TABLE1_hg38): drops known RNA editing sites
(same chrom, position, ref and edited base).
RNA-seq artifact filters: indels longer than --max_indel_len (usually unannotated splice
junctions read as deletions) and clusters of calls within --cluster_dist bp (editing
clusters, misalignment) are dropped.
"""
import argparse
import gzip
import sys


def sample_values(fmt, col):
    return dict(zip(fmt.split(':'), col.split(':')))


def nums(s):
    return [float(x) for x in s.split(',')]


def keep(f, tum, nor, a):
    if f[6] != 'PASS':
        return False
    info = dict(x.split('=', 1) for x in f[7].split(';') if '=' in x)
    # POPAF = -log10(population AF); AF < max_pop_af  <=>  POPAF > -log10(max_pop_af)
    if 'POPAF' in info and max(10 ** -p for p in nums(info['POPAF'])) >= a.max_pop_af:
        return False
    if max(abs(len(f[3]) - len(alt)) for alt in f[4].split(',')) > a.max_indel_len:
        return False
    t_alt = max(nums(tum['AD'])[1:])
    n_alt = max(nums(nor['AD'])[1:])
    return (float(tum['DP']) >= a.min_tumor_dp
            and float(nor['DP']) >= a.min_normal_dp
            and t_alt >= a.min_tumor_alt
            and max(nums(tum['AF'])) >= a.min_tumor_af
            and (max(nums(nor['AF'])) <= a.max_normal_af or n_alt <= a.max_normal_alt))


def editing_hits(path, wanted):
    """Return the (chrom, pos, ref, alt) keys from `wanted` that are known REDIportal editing sites.
    REDIportal columns: Accession, Region, Position, Ref, Ed (genomic strand)."""
    hits = set()
    with gzip.open(path, 'rt') as fh:
        next(fh)  # header
        for line in fh:
            f = line.split('\t', 5)
            key = (f[1], f[2], f[3], f[4])
            if key in wanted:
                hits.add(key)
    return hits


def clustered(keys, dist):
    """Keys (chrom, pos, ...) of calls that have another call within `dist` bp on the same chrom."""
    hits = set()
    ordered = sorted(keys, key=lambda k: (k[0], int(k[1])))
    for prev, cur in zip(ordered, ordered[1:]):
        if prev[0] == cur[0] and int(cur[1]) - int(prev[1]) <= dist:
            hits.update((prev, cur))
    return hits


def demo():
    k = [('chr1', '100', 'A', 'G'), ('chr1', '105', 'T', 'C'), ('chr1', '200', 'C', 'T'), ('chr2', '103', 'G', 'A')]
    assert clustered(k, 10) == {k[0], k[1]}
    assert clustered(k, 4) == set()
    print('pass_filter demo OK')


def main(argv=None):
    p = argparse.ArgumentParser()
    if '--demo' in (argv or sys.argv[1:]):
        return demo()
    p.add_argument('vcf')
    p.add_argument('--tumor', required=True)
    p.add_argument('--normal', required=True)
    p.add_argument('--min_tumor_dp', type=float, default=20)
    p.add_argument('--min_normal_dp', type=float, default=20)
    p.add_argument('--min_tumor_alt', type=float, default=5)
    p.add_argument('--min_tumor_af', type=float, default=0.05)
    p.add_argument('--max_normal_af', type=float, default=0.02)
    p.add_argument('--max_normal_alt', type=float, default=1)
    p.add_argument('--max_pop_af', type=float, default=0.001)
    p.add_argument('--max_indel_len', type=int, default=20)
    p.add_argument('--cluster_dist', type=int, default=10, help='0 disables the cluster filter')
    p.add_argument('--editing_sites', help='REDIportal TABLE1_hg38*.txt.gz (optional)')
    a = p.parse_args(argv)

    kept = total = 0
    out = []  # (key, line); written at the end so editing sites can be dropped
    with gzip.open(a.vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                sys.stdout.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                sys.stdout.write(f'##pass_filter={" ".join(sys.argv[1:])}\n' + line)
                ti, ni = f.index(a.tumor), f.index(a.normal)
                continue
            total += 1
            try:
                ok = keep(f, sample_values(f[8], f[ti]), sample_values(f[8], f[ni]), a)
            except (KeyError, ValueError):  # missing field or '.'
                ok = False
            if ok:
                kept += 1
                out.append(((f[0], f[1], f[3], f[4]), line))
    keys = {k for k, _ in out}
    edits = editing_hits(a.editing_sites, keys) if a.editing_sites else set()
    clusters = clustered(keys, a.cluster_dist) - edits if a.cluster_dist else set()
    written = 0
    for k, line in out:
        if k not in edits and k not in clusters:
            sys.stdout.write(line)
            written += 1
    print(f'{a.vcf}: {kept} of {total} variants pass cutoffs; {len(edits)} removed as known RNA editing sites; '
          f'{len(clusters)} removed as clustered (within {a.cluster_dist} bp); {written} written', file=sys.stderr)


if __name__ == '__main__':
    main()
