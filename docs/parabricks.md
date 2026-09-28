# Optional Parabricks Mutect2 backend

This fork of [genepi/mtdna-server-2](https://github.com/genepi/mtdna-server-2)
adds an opt-in raw Mutect2 caller using NVIDIA Parabricks **4.7.1-1**.
The default `--mutect2_backend gatk` retains the GATK calling command.
Upstream base: `7e41e355c4fd6ebac7f860b5011e55b1c2ce0569` (v2.1.16).

```bash
nextflow run phonegor95/mtdna-server-2 -r main \
  -profile singularity,parabricks \
  --project example --files '/data/*.bam' --output results \
  --mode fusion --detection_limit 0.1 --parabricks_gpus 1
```

Pin `-r` to the desired fork commit for reproducible deployments. A working
NVIDIA driver and GPU-aware Singularity/Apptainer or Docker installation are
required. `-profile docker,parabricks` uses Docker's `--gpus` option.
The default caller requests 8 CPU threads and 48 GB host memory; override
`parabricks_cpus` and `parabricks_memory` for the deployment. The raw caller
is limited to one concurrent task by default. On a shared host, scope
`CUDA_VISIBLE_DEVICES` to reserved devices or provide a site launcher;
Nextflow local execution alone does not reserve devices across workflows.

`--parabricks_container /path/to/parabricks.sif` selects a prepared image.
Other stages continue to use the upstream supporting image. GATK
`FilterMutectCalls`, normalization, the allele-fraction cutoff, MUTSERVE,
fusion merging, annotation and report generation retain their existing commands.
Fusion selects Mutect2 indels and MUTSERVE SNVs as upstream does.
`--mode mutserve` does not invoke either Mutect2 backend.

## Compatibility

The adapter reads the BAM's SM tags and requires exactly one sample across
all read groups. It forwards the mitochondrial contig, `baseQ` and disabled
read-start downsampling, and requires the raw VCF, index and `.stats` output
before GATK filtering. It does **not** enable `--mitochondria-mode`, because
the upstream GATK invocation does not use that option.

Parabricks does not expose upstream's `--callable-depth 6` setting. Its
callable-site statistics and consequently filtering may differ at low depth.
Treat the backend as experimental; synthetic tests do not establish clinical
or patient-level equivalence or a speedup. See NVIDIA's
[mutectcaller documentation](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/mutectcaller).

## GenDecoder integration

GenDecoder's `gpu.mtdna_mutect2: true` selects this fork at its pinned commit.
Its host launcher leases GPUs from the same pool as Sarek and QC tasks.
The fork remains usable without GenDecoder. For other host launchers,
`parabricks_runner` accepts a list of argv strings prepended to `pbrun`;
`parabricks_python` and `parabricks_samtools` select host tools. Override the
`PARABRICKS_MUTECT2` process's `container` to `null` and `containerOptions`
to an empty string when the launcher manages the container itself.
`parabricks_signature` can include a runtime/configuration hash for Nextflow's cache.

## Validation

```bash
python3 -m unittest discover -s tests/parabricks -p 'test_*.py'
python3 tests/parabricks/validate.py \
  --outdir /scratch/mt-validation \
  --gatk-image /images/gatk-with-samtools-bcftools-tabix.sif \
  --parabricks-image /images/parabricks.sif
```

The validation script requires `pysam`, Nextflow 25.10.4, samtools and
Singularity. The supporting GATK image must also contain bcftools and tabix.
It generates synthetic reads without accessing patient data, runs both real
callers and GATK filtering, checks sample identity, calls, allele fractions
and filter status, and exercises the fusion indel-selection module.
Outputs include traces, logs and a `validation.json` receipt. Full fusion
reports and patient-level concordance are separate validation steps.

Validated on 2026-09-28 with Parabricks 4.7.1-1, GATK 4.7.0.0,
Nextflow 25.10.4 and an RTX PRO 5000 72GB Blackwell GPU. Both native
Singularity execution and a host GPU launcher passed. Synthetic SNVs at
5% and 9% were excluded by the 10% cutoff; 11%, 20%, 50% and 100% sites
were retained with matching allele fractions and filter status. A 9%
insertion was excluded; 20% insertion/deletion calls passed and were
retained by fusion selection. The complete workflow also passed a
Nextflow preview; this is a compilation check, not a full report run.
