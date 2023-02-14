```
module load singularity/3.1.1 nextflow/21.10.6

nextflow run variant-analysis.nf -profile beartooth

nextflow -bg run variant-analysis.nf -progile beartooth

nextflow -bg run variant-analysis.nf -progile beartooth --max_snp_missingness 18

nextflow -bg run variant-analysis.nf -progile beartooth --max_snp_missingness 9

nextflow run variant-analysis.nf -profile beartooth && nextflow clean -f -q

```
