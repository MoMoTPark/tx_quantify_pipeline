### Transcript/gene region/exon RNA-seq data quantitation pipeline

#### Pipeline description:
This pipeline is prepared to quantitate custom designed experiments that contain pair-end sequencing data with overlapping paired-end reads. The overlapping pair-end data are merged into one long read and treated as ground truth representation of the covering regions. Then this data is used for quantifying (with normalisation) a select region of the transcriptome.

#### Relavent code:
`snakefile` contains the pipeline in Snakemake format.