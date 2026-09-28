#!/usr/bin/env python3
"""Synthetic caller/filter comparison; requires pysam, Nextflow, Singularity and a GPU."""
import argparse
import csv
import json
import os
from pathlib import Path
import random
import subprocess
import pysam

ROOT = Path(__file__).resolve().parents[2]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--gatk-image',required=True)
    p.add_argument('--parabricks-image',required=True)
    p.add_argument('--nextflow',default='nextflow')
    p.add_argument('--samtools',default='samtools')
    p.add_argument('--runner-json',default='[]')
    a=p.parse_args()
    out=a.outdir.resolve();out.mkdir(parents=True,exist_ok=False)
    rng=random.Random(178)
    seq=''.join(rng.choices('ACGT',k=16569))
    # Avoid a repeat around the deletion so its normalized position is fixed.
    seq=seq[:2998]+'ACGT'+seq[3002:]
    ref=out/'ref.fasta';ref.write_text('>chrM\n'+seq+'\n')
    pysam.faidx(str(ref))
    with (out/'ref.dict').open('w') as f:subprocess.run([a.samtools,'dict',str(ref)],stdout=f,check=True)
    sites={1000:0.20,2000:0.09,3000:0.20,5000:0.05,7000:0.09,9000:0.11,11000:0.20,13000:0.50,15000:1.0}
    expected={}
    with pysam.AlignmentFile(str(out/'unsorted.bam'),'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chrM','LN':len(seq)}],'RG':[{'ID':'rg1','SM':'synthetic_mt','PL':'ILLUMINA','LB':'synthetic'}]}) as bam:
        for pos,fraction in sites.items():
            alt=next(b for b in 'ACGT' if b!=seq[pos-1])
            kind='ins' if pos in (1000,2000) else ('del' if pos==3000 else 'snv')
            expected[pos]=dict(fraction=fraction,ref=seq[pos-1:pos+1] if kind=='del' else seq[pos-1],
                               alt=seq[pos-1] if kind=='del' else (seq[pos-1]+alt if kind=='ins' else alt),kind=kind)
            choices=[True]*int(1000*fraction)+[False]*(1000-int(1000*fraction));rng.shuffle(choices)
            for n,is_alt in enumerate(choices):
                start=pos-1-rng.randrange(35,115);bases=list(seq[start:start+150])
                cigar='150M'
                if is_alt:
                    offset=pos-start
                    if kind=='ins':
                        bases=(bases[:offset]+[alt]+bases[offset:])[:150]
                        cigar=f'{offset}M1I{149-offset}M'
                    elif kind=='del':
                        bases=list(seq[start:start+151]);del bases[offset]
                        cigar=f'{offset}M1D{150-offset}M'
                    else:bases[pos-1-start]=alt
                read=pysam.AlignedSegment();read.query_name=f'read_{pos}_{n}';read.query_sequence=''.join(bases)
                read.flag=16 if n%2 else 0;read.reference_id=0;read.reference_start=start
                read.mapping_quality=60;read.cigarstring=cigar;read.query_qualities=pysam.qualitystring_to_array('I'*150);read.set_tag('RG','rg1');bam.write(read)
    bam=out/'synthetic_mt.bam';pysam.sort('-o',str(bam),str(out/'unsorted.bam'));pysam.index(str(bam));(out/'unsorted.bam').unlink()
    results={}
    for backend in ('gatk','parabricks'):
        launch=out/backend;launch.mkdir()
        params=dict(mutect2_backend=backend,baseQ=20,detection_limit=0.1,
                    parabricks_container=a.parabricks_image,parabricks_gpus=1,parabricks_cpus=8,
                    test_bam=str(bam),test_reference=str(ref),test_fai=str(ref)+'.fai',test_dict=str(out/'ref.dict'),
                    parabricks_runner=json.loads(a.runner_json),parabricks_samtools=a.samtools if json.loads(a.runner_json) else 'samtools')
        (launch/'params.json').write_text(json.dumps(params,indent=2))
        config="""process.executor = 'local'
process.container = %s
process.memory = '4 GB'
process.cpus = 2
process.errorStrategy = 'terminate'
process.maxRetries = 0
docker.enabled = false
singularity.enabled = true
singularity.autoMounts = true
tower.enabled = false
notification.enabled = false
process { withName: FILTER_VARIANTS { publishDir = [path: "%s", mode: 'copy', pattern: '*.filtered.txt'] } }
process { withName: MUTECT2 { publishDir = [path: "%s", mode: 'copy', pattern: 'synthetic_mt.bam.vcf.gz*'] } }
"""%(json.dumps(a.gatk_image),str(launch/'results'),str(launch/'results'))
        # bam_file.baseName for synthetic_mt.bam is synthetic_mt (no .bam).
        config=config.replace('synthetic_mt.bam.vcf.gz*','synthetic_mt.vcf.gz*')
        if a.runner_json!='[]':config+="process { withName: PARABRICKS_MUTECT2 { container = null; containerOptions = '' } }\n"
        (launch/'test.config').write_text(config)
        with (launch/'nextflow.log').open('w') as log:
            subprocess.run([a.nextflow,'run',str(ROOT/'tests/parabricks/main.nf'),'-c',str(ROOT/'nextflow.config'),'-c',str(launch/'test.config'),'-params-file',str(launch/'params.json'),'-work-dir',str(launch/'work'),'-with-trace',str(launch/'trace.tsv')],cwd=launch,env={**os.environ,'NXF_VER':'25.10.4','NXF_OFFLINE':'true'},stdout=log,stderr=subprocess.STDOUT,check=True)
        with pysam.VariantFile(str(launch/'results/synthetic_mt.vcf.gz')) as vcf:
            assert list(vcf.header.samples)==['synthetic_mt']
            results[backend]={r.pos:dict(ref=r.ref,alt=list(r.alts),af=r.samples['synthetic_mt']['AF'][0],filters=list(r.filter)) for r in vcf}
        fusion_file=next((launch/'results').glob('*.filtered.txt'))
        with fusion_file.open() as f:
            fusion_rows=list(csv.DictReader(f,delimiter='\t'))
        assert {int(row['Pos']) for row in fusion_rows}=={1000,3000}, fusion_rows
        assert set(results[backend])=={pos for pos,v in expected.items() if v['fraction']>=0.1}, results[backend]
        for pos,r in results[backend].items():
            assert r['ref']==expected[pos]['ref'] and r['alt']==[expected[pos]['alt']]
            assert abs(r['af']-expected[pos]['fraction'])<0.025
    for pos in results['gatk']:
        assert results['gatk'][pos]['filters']==results['parabricks'][pos]['filters']
        assert abs(results['gatk'][pos]['af']-results['parabricks'][pos]['af'])<0.005
    (out/'validation.json').write_text(json.dumps(dict(status='passed',scope='synthetic SNVs/indels, GATK filtering and fusion indel selection; not full pipeline or patient validation',cutoff=0.1,results=results),indent=2)+'\n')
    print((out/'validation.json').read_text())
if __name__=='__main__':main()
