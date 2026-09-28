#!/usr/bin/env python3
"""Run the optional raw caller; GATK filtering stays in the Nextflow module."""
import argparse
import json
from pathlib import Path
import subprocess


def sample_name(header):
    samples = set()
    for line in header.splitlines():
        if line.startswith('@RG\t'):
            names = [field[3:] for field in line.split('\t') if field.startswith('SM:')]
            if len(names) != 1 or not names[0]:
                raise ValueError('Every read group must contain one nonempty SM tag')
            samples.update(names)
    if len(samples) != 1:
        raise ValueError('Parabricks requires exactly one sample across all read groups')
    return samples.pop()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bam', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--contig', required=True)
    p.add_argument('--base-quality', type=int, required=True)
    p.add_argument('--gpus', type=int, default=1)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--runner-json', default='[]')
    a = p.parse_args()
    runner = json.loads(a.runner_json)
    if not isinstance(runner, list) or any(not isinstance(x, str) for x in runner):
        p.error('--runner-json must contain an argv list')
    if a.gpus < 1 or a.threads < 1:
        p.error('GPU and thread counts must be positive')
    header = subprocess.check_output([a.samtools, 'view', '-H', a.bam], text=True)
    sample = sample_name(header)
    subprocess.run([a.samtools, 'index', a.bam], check=True)
    command = ['pbrun', 'mutectcaller', '--ref', a.reference,
               '--in-tumor-bam', a.bam, '--tumor-name', sample,
               '--out-vcf', 'raw.vcf.gz', '--interval', a.contig,
               '--min-base-quality-score', str(a.base_quality),
               '--mutectcaller-options=-max-reads-per-alignment-start 0',
               '--num-gpus', str(a.gpus),
               '--num-htvc-threads', str(max(1, min(5, a.threads // a.gpus))),
               '--preserve-file-symlinks']
    # Upstream does not use --mitochondria-mode; do not change calling thresholds.
    subprocess.run(runner + command, check=True)
    for name in ('raw.vcf.gz', 'raw.vcf.gz.stats', 'raw.vcf.gz.tbi'):
        if not Path(name).is_file() or Path(name).stat().st_size == 0:
            raise RuntimeError('Missing required caller output: ' + name)


if __name__ == '__main__':
    main()
