```

module load arcc/1.0 gcc/12.2.0 git/2.37.0 nextflow/22.10.4 singularity/3.10.3
NXF_SINGULARITY_CACHEDIR=$(pwd)/.singularity/cache

nextflow run plugnplay-popgen.nf -c config/beartooth.config

nextflow -bg run variant-analysis.nf -progile beartooth

nextflow -bg run variant-analysis.nf -progile beartooth --max_snp_missingness 18

nextflow -bg run variant-analysis.nf -progile beartooth --max_snp_missingness 9

nextflow run variant-analysis.nf -profile beartooth && nextflow clean -f -q

```
