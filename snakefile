## Merge and process overlapping paired-end short read RNA-seq data to quantify a section of given gene
## This pipeline generates the outputs required for further downstream quantitation using R, etc.

# Note this snakefile should be executed from a directory (e.g., code) within root analysis directory.
# Hence, relative path is provided with '../' prefix 

IDS, = glob_wildcards("../data/{id}_R1.fastq.gz")

rule all:
    input:
        expand("../fastq_files/{id}_dedup_R1.fastq.gz", id=IDS),
        expand("../fastq_files/{id}_dedup_R2.fastq.gz", id=IDS),
        expand("../fastq_files/{id}_dedup_R1_header.fastq.gz", id=IDS),
        expand("../fastq_files/{id}_dedup_R2_header.fastq.gz", id=IDS),
        expand("../merged_reads/{id}_merged_reads.fastq.gz", id=IDS),
        expand("../star_align2/{id}.Aligned.out.bam", id=IDS),
        expand("../salmon_quant/{id}/quant.sf", id=IDS)

# Generate QC html from raw reads for agent trim
# -A: Disable adaptor trimming
# -G: Disable trim poly-G
# -Q: Disable quality read filtering
# -L: Disable length filtering
# -w: number of threads
# -h: html output
rule fastp_qc:
    input:
        read_1 = "../data/{id}_R1.fastq.gz",
        read_2 = "../data/{id}_R2.fastq.gz"
    output:
        html = "../fastp_qc/{id}_fastp_qc.html",
        json = "../fastp_qc/{id}_fastp_qc.json"
    log: "../logs/{id}_fastp_qc.log"
    benchmark: "../benchmarks/{id}_fastp_qc.benchmark"
    shell:
        '''fastp -i {input.read_1} -I {input.read_2} -A -G -Q -L -w 4 -h {output.html} -j {output.json} 2> {log}'''

# -v2 specifies HS2 library prep protocol
rule agent_trim:
    input:
        read_1 = "../data/{id}_R1.fastq.gz",
        read_2 = "../data/{id}_R2.fastq.gz",
        html = "../fastp_qc/{id}_fastp_qc.html"
    output:
        agent_out = "../agent_trim/{id}_MBC.txt.gz",
        read_1 = "../agent_trim/{id}_R1.fastq.gz",
        read_2 = "../agent_trim/{id}_R2.fastq.gz"
    params:
        out_prefix = "../agent_trim/{id}"
    log: "../logs/{id}_agent_trim.log"
    benchmark: "../benchmarks/{id}_agent_trim.benchmark"
    shell:
        '''agent trim -v2 -fq1 {input.read_1} -fq2 {input.read_2} -polyG 3 -qualityTrimming 15 -out {params.out_prefix} 2> {log}'''

# Trim reads with fastp with the following parameters
# -q, the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
# -u, how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
# -y, enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
# -Y, the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
rule fastp_trim:
    input:
        read_1 = "../agent_trim/{id}_R1.fastq.gz",
        read_2 = "../agent_trim/{id}_R2.fastq.gz"
    output:
        html = "../fastp/{id}_fastp.html",
        json = "../fastp/{id}_fastp.json",
        trimmed_read_1 = "../trimmed_reads/{id}_R1_qctrim.fastq.gz",
        trimmed_read_2 = "../trimmed_reads/{id}_R2_qctrim.fastq.gz"
    log: "../logs/{id}_fastp_trim.log"
    benchmark: "../benchmarks/{id}_fastp_trim.benchmark"
	shell:
		'''fastp -i {input.read_1} -I {input.read_2} -q 15 -u 40 -y -Y 30 -h {output.html} -j {output.json} -o {output.trimmed_read_1} -O {output.trimmed_read_2} 2> {log}'''

# Preserve fastq headers with barcode information included and output as bam with BC tags included 
rule samtools_import:
    input:
        read_1 = "../trimmed_reads/{id}_R1_qctrim.fastq.gz",
        read_2 = "../trimmed_reads/{id}_R2_qctrim.fastq.gz"
    output:
        bam = "../bam_files/{id}.bam"
    log: "../logs/{id}_samtools_import.log"
    benchmark: "../benchmarks/{id}_samtools_import.benchmark"
    shell:
        '''samtools import -i -1 {input.read_1} -2 {input.read_2} -o {output.bam} -O BAM -@ 4 -T '*' 2> {log}'''

# First alignment with STAR
rule star_align_1:
    input:
        bam_reads = "../bam_files/{id}.bam",
        genome_index = '../star_index/Log.out'
    output:
        star_out = '../star_align1/{id}.SJ.out.tab',
        alignment = '../star_align1/{id}.Aligned.out.bam'
    params:
        genome_index = '../star_index',
        outfile_prefix = '../star_align1/{id}.',
        outfile_tmp_dir = '../tmp/{id}'
    log: "../logs/{id}_star_align1.log"
    benchmark: "../benchmarks/{id}_star_align1.benchmark"
    resources:
        load = 25 
    shell:
        '''	STAR \
        --runThreadN 8 \
        --genomeDir {params.genome_index} \
        --readFilesIn {input.bam_reads} \
        --readFilesType SAM PE \
        --readFilesCommand samtools view \
        --outFilterType BySJout \
        --outFileNamePrefix {params.outfile_prefix} \
        --outSAMtype BAM Unsorted \
        2> {log}'''

# Sort BAM files 
rule sort_bam:
	input: '../star_align1/{id}.Aligned.out.bam'
	output: '../star_align1/{id}.sorted.bam'
	log: "../logs/{id}_pos_sort_bam.log"
	benchmark: "../benchmarks/{id}_pos_sort_bam.benchmark"
	shell:
		'''samtools sort -@ 4 -o {output} {input} 2> {log}'''	

# Query sort aligned bam file
rule qsort_bam:
    input: "../star_align1/{id}.sorted.bam"
    output: "../star_align1/{id}.qsorted.bam"
    log: "../logs/{id}_qsort_bam.log"
    benchmark: "../benchmarks/{id}_qsort_bam.benchmark"
    shell:
        '''samtools sort -@ 4 -n -o {output} -O bam {input} 2> {log}'''

# Sort MBC output from agent trim and generate fastq.gz output
rule fgbio_sort_fastq:
    input:
        jar_file = "fgbio/share/fgbio/fgbio.jar",
        mbc = "../agent_trim/{id}_MBC.txt.gz"
    output:
        mbc_fastq = "../dedup/{id}_MBC.fastq.gz"
    log: "../logs/{id}_fgbio_sort_fastq.log"
    benchmark: "../benchmarks/{id}_fgbio_sort_fastq.benchmark"
    shell:
        '''java -jar {input.jar_file} SortFastq -i {input.mbc} -o {output.mbc_fastq} 2> {log}'''

# Added a step to change the reference name to something unique as locatit does not like the naming of the files
# Change reference name from '3:193615674-193631665' to 'hola'
rule change_bam_name_1:
    input:
        queryname_sorted_bam = "../star_align1/{id}.qsorted.bam"
    output:
        bam = "../star_align1/{id}.qsorted_locatit.bam"
    log: "../logs/{id}_change_bam_name_1.log"
    benchmark: "../benchmarks/{id}_change_bam_name_1.benchmark"
    shell:
        """
        samtools view -h {input.queryname_sorted_bam} | perl -pe 's/3\:193615674\-193631665/hola/g' | samtools view -Sb - > {output.bam}
        """

# De-duplicate reads and generate deduplicated bam files
# -i: Switch to turn on SureSelect and turn off HaloPlex
# -R: Remove duplicates(SAM flag 0x400), filtered reads(SAM flag 0x200), secondary(SAM flag 0x100) and supplementary(SAM flag 0x800) reads
# -C: For SureSelect it should be always on. If -i is specified it sets -C internally.
# -q: Reads having barcodes with base quality less than specified threshold will be filtered. Default is 25. Range: 0-45
# -L: Allows LocatIt to discard non-matching index2 reads as it processes the input file instead of buffering them. This saves a lot of memory.
# -m: Minimum number of read pairs associated with a barcode. Default is 1. Range: >=1
# -X: Location of temporary files used to store overflow of matches. Temporary files will be deleted at program exit.

rule agent_locatit:
    input:
        queryname_sorted_bam = "../star_align1/{id}.qsorted_locatit.bam",
        mbc_fastq = "../dedup/{id}_MBC.fastq.gz"
    output: "../dedup/{id}_dedup_locatit.bam"
    params:
        tmp = "../tmp/{id}"
    log: "../logs/{id}_agent_locatit.log"
    benchmark: "../benchmarks/{id}_agent_locatit.benchmark"
    shell:
        '''agent -Xmx20G locatit -i -c 2500 -C -q 0 -U -L -m 1 -d 1 -R \
        -X {params.tmp} -PM:xm,Q:xq,q:nQ,r:nR \
        -IB -OB -o {output} {input.queryname_sorted_bam} {input.mbc_fastq} 2> {log}'''

# Added a step to change the reference name back to 3:193615674-193631665
rule change_bam_name_2:
    input: "../dedup/{id}_dedup_locatit.bam"
    output: "../dedup/{id}_dedup.bam"
    log: "../logs/{id}_change_bam_name_2.log"
    benchmark: "../benchmarks/{id}_change_bam_name_1.benchmark"
    shell:
        """
        samtools view -h {input} | perl -pe 's/hola/3\:193615674\-193631665/g' | samtools view -Sb - > {output}
        """

# Sort bam files
rule sort_dedup_bam:
	input: "../dedup/{id}_dedup.bam"
	output: "../dedup/{id}_dedup_sorted.bam"
	log: "../logs/{id}_dedup_pos_sort_bam.log"
	benchmark: "../benchmarks/{id}_dedup_pos_sort_bam.benchmark"
	shell:
		'''samtools sort -@ 4 -o {output} {input} 2> {log}'''

# Query sort bam files
rule qsort_dedup_bam:
    input: "../dedup/{id}_dedup.bam"
    output: "../dedup/{id}_dedup.qsorted.bam"
    log: "../logs/{id}_qsort_dedup_bam.log"
    benchmark: "../benchmarks/{id}_qsort_dedup_bam.benchmark"
    shell:
        '''samtools sort -m 8G -@ 4 -n -o {output} {input} 2> {log}'''

# Convert BAM to FASTQ
rule bam_to_fastq:
    input:
        bam = "../dedup/{id}_dedup.qsorted.bam"
    output:
        read_1 = "../fastq_files/{id}_dedup_R1.fastq.gz",
        read_2 = "../fastq_files/{id}_dedup_R2.fastq.gz"
    log: "../logs/{id}_bam2fastq.log"
	benchmark: "../benchmarks/{id}_bam2fastq.benchmark"
	shell:
		'''
        bedtools bamtofastq -i {input.bam} -fq {output.read_1} -fq2 {output.read_2} 2> {log}
        '''

# Remove /1 and /2 from paired fastq header for seqkit pair compatibility
rule modify_fastq_headers:
    input:
        read_1 = "../fastq_files/{id}_dedup_R1.fastq.gz",
        read_2 = "../fastq_files/{id}_dedup_R2.fastq.gz"
    output:
        header_1 = "../fastq_files/{id}_dedup_R1_header.fastq.gz",
        header_2 = "../fastq_files/{id}_dedup_R2_header.fastq.gz"
    log: "../logs/{id}_modify_fastq_headers.log"
    benchmark: "../benchmarks/{id}_modify_fastq_headers.benchmark"
    shell:
        """
        sed '/^@.*\\/1$/ s/\\/1$//' {input.read_1} > {output.header_1};
        sed '/^@.*\\/2$/ s/\\/2$//' {input.read_2} > {output.header_2} 2> {log}
        """

# Check fastq file integrity
rule seqkit_sana_fastq_check:
    input:
        read_1 = "../fastq_files/{id}_dedup_R1_header.fastq.gz",
        read_2 = "../fastq_files/{id}_dedup_R2_header.fastq.gz"
    output:
        checked_1 = "../fastq_files/{id}_dedup_R1_checked.fastq.gz",
        checked_2 = "../fastq_files/{id}_dedup_R2_checked.fastq.gz"
    log: "../logs/{id}_seqkit_sana_fastq_check.log"
    benchmark: "../benchmarks/{id}_seqkit_sana_fastq_check.benchmark"
    shell:
        """
        ~/tools/seqkit sana {input.read_1} -o {output.checked_1} 2> {log};
        ~/tools/seqkit sana {input.read_2} -o {output.checked_2}
        """

# Remove reads with no mate pair
rule remove_unpaired_reads:
    input:
        read_1 = "../fastq_files/{id}_dedup_R1_checked.fastq.gz",
        read_2 = "../fastq_files/{id}_dedup_R2_checked.fastq.gz"
    output:
        paired_1 = "../fastq_files/{id}_dedup_R1_checked.paired.fastq.gz",
        paired_2 = "../fastq_files/{id}_dedup_R2_checked.paired.fastq.gz"
    log: "../logs/{id}_remove_unpaired_reads.log"
    benchmark: "../benchmarks/{id}_remove_unpaired_reads.benchmark"
    shell:
        """
        ~/tools/seqkit pair -1 {input.read_1} -2 {input.read_2} 2> {log}
        """

# Merge overlapping paired-end reads with NGmerge
rule ngmerge:
    input:
        read_1 = "../fastq_files/{id}_dedup_R1_checked.paired.fastq.gz",
        read_2 = "../fastq_files/{id}_dedup_R2_checked.paired.fastq.gz"
    output:
        merge_log = "../merged_reads/{id}_merge.log",
        merge_alignment_log = "../merged_reads/{id}_merge_alignment.log",
        failed_reads = "../merged_reads/{id}_failed_1.fastq.gz",
        merged_reads = "../merged_reads/{id}_merged_reads.fastq.gz"
    params:
        failed = "../merged_reads/{id}_failed"
    log: "../logs/{id}_ngmerge.log"
    benchmark: "../benchmarks/{id}_ngmerg.benchmark"
    shell:
        '''NGmerge -1 {input.read_1} -2 {input.read_2} -o {output.merged_reads} -l {output.merge_log} -j {output.merge_alignment_log} -f {params.failed} -n 4 -z 2> {log}'''

# Merged read alignment with STAR
rule star_align_2:
    input:
        merged_reads = "../merged_reads/{id}_merged_reads.fastq.gz",
        genome_index = '../star_index/Log.out'
    output:
        star_out = '../star_align2/{id}.SJ.out.tab',
        alignment = '../star_align2/{id}.Aligned.out.bam'
    params:
        genome_index = '../star_index',
        outfile_prefix = '../star_align2/{id}.',
        outfile_tmp_dir = '../tmp/{id}'
    log: "../logs/{id}_star_align2.log"
    benchmark: "../benchmarks/{id}_star_align2.benchmark"
    resources:
        load = 25 
    shell:
        '''	STAR \
        --runThreadN 8 \
        --genomeDir {params.genome_index} \
        --readFilesIn {input.merged_reads} \
        --readFilesType Fastx \
        --readFilesCommand gunzip -c \
        --outFilterType BySJout \
        --readQualityScoreBase 28 \
        --outFileNamePrefix {params.outfile_prefix} \
        --quantMode TranscriptomeSAM GeneCounts \
        --outBAMcompression 9 \
        --outSAMtype BAM Unsorted \
        2> {log}'''

### The following two steps only needed if we suspect there are unique splice junctions in our dataset
# Prepare for realigning reads to assembled contigs
# rule splice_junction_filter:
#     input:
#         splice_junctions = '../star_align2/{id}.SJ.out.tab'
#     output:
#         splice_junctions_filtered = '../star_align2/{id}.SJ.out.tab.filter'
#     log: "../logs/{id}_splice_junction_filter.log"
#     benchmark: "../benchmarks/{id}_splice_junction_filter.benchmark"
#     shell:
#         '''awk '$1~/chr[1-2XY]/ && $6==0 && $5>0 && $7>0' {input.splice_junctions} > {output.splice_junctions_filtered} 2> {log}'''

# Align reads to splice junctions and reference  
# rule star_align_3:
#     input:
#         merged_reads = "{id}_merged_reads.fastq",
#         genome_index = '../star_index/Log.out',
#         splice_junctions_filtered = '../star_align2/{id}.SJ.out.tab.filter'
#     output:
#         star_aligned3 = '../star_align3/{id}.Aligned.out.bam'
#     params:
#         genome_index = '../star_index',
#         outfile_prefix = '../star_align3/{id}.'
#     log: "../logs/{id}_star_align3.log"
#     benchmark: "../benchmarks/{id}_star_align3.benchmark"
#     resources:
#         load = 33 
#     shell:
#         '''STAR \
#         --runThreadN 8 \
#         --genomeDir {params.genome_index} \
# 	    --readFilesIn {input.merged_reads} \
#         --readFilesType Fastx \
# 	    --outFilterType BySJout \
#         --outFilterMultimapNmax 20 \
#         --alignSJoverhangMin 8 \
#         --alignSJDBoverhangMin 1 \
#         --outFilterMismatchNmax 999 \
#         --outFilterMismatchNoverReadLmax 0.04 \
#         --alignIntronMin 20 \
#         --alignIntronMax 1000000 \
#         --alignMatesGapMax 1000000 \
# 	    --sjdbFileChrStartEnd {input.splice_junctions_filtered} \
# 	    --quantMode TranscriptomeSAM GeneCounts \
#         --limitSjdbInsertNsj 20000000 \
# 	    --outBAMcompression 9 \
# 	    --outFileNamePrefix {params.outfile_prefix} \
# 	    --outSAMtype BAM Unsorted	\
#         2> {log}'''

        # ENCODE alignment parameters
        # --outFilterMultimapNmax 20 \
        # --alignSJoverhangMin 8 \
        # --alignSJDBoverhangMin 1 \
        # --outFilterMismatchNmax 999 \
        # --outFilterMismatchNoverReadLmax 0.04 \
        # --alignIntronMin 20 \
        # --alignIntronMax 1000000 \
        # --alignMatesGapMax 1000000 \

# Sort BAM files 
rule sort_bam_2:
	input: '../star_align2/{id}.Aligned.out.bam'
	output: '../star_align2/{id}.sorted.bam'
	log: "../logs/{id}_pos_sort_bam_2.log"
	benchmark: "../benchmarks/{id}_pos_sort_bam_2.benchmark"
	shell:
		'''samtools sort -@ 4 -o {output} {input} 2> {log}'''	

# Query sort aligned bam file
rule qsort_bam_2:
    input: "../star_align2/{id}.sorted.bam"
    output: "../star_align2/{id}.qsorted.bam"
    log: "../logs/{id}_qsort_bam_2.log"
    benchmark: "../benchmarks/{id}_qsort_bam_2.benchmark"
    shell:
        '''samtools sort -@ 4 -n -o {output} -O bam {input} 2> {log}'''

# Convert bam to single read fasq files
rule bam_to_fastq_2:
    input:
        bam = "../star_align2/{id}.qsorted.bam"
    output:
        reads = "../fastq_files/{id}_salmon.fastq.gz"
    log: "../logs/{id}_bam2fastq_2.log"
    benchmark: "../benchmarks/{id}_bam2fastq_2.benchmark"
    shell:
        """
        bedtools bamtofastq -i {input.bam} -fq {output.reads} \
        2> {log}
        """

# Quantify with Salmon
rule salmon_quant:
    input:
        reads = "../fastq_files/{id}_salmon.fastq.gz",
        index = "../salmon_index/"
    output:
        quant = "../salmon_quant/{id}/quant.sf"
    params:
        outfolder = "../salmon_quant/{id}"
    log: "../logs/{id}_salmon_quant.log"
    benchmark: "../benchmarks/{id}_salmon_quant.benchmark"
    shell:
        """
        salmon quant -i {input.index} -l A -r {input.reads} --validateMappings --threads 4 -o {params.outfolder} \
        2> {log}
        """

# Index BAM files
# rule index_bam:
# 	input: "../star_align2/{id}.sorted.bam"
# 	output: "../star_align2/{id}.sorted.bam.bai"
# 	log: "../logs/{id}_index_bam.log"
# 	benchmark: "../benchmarks/{id}_index_bam.benchmark"
# 	shell:
# 		'''samtools index -@ 4 {input} 2> {log}'''	

# # Run QoRTs
# rule run_qorts:
# 	input:
# 		bam = "../star_align2/{id}.sorted.bam",
# 		gtf = "../ref/OPA1_3_193615674_193631665.gtf"
# 	output:
# 		qorts = '../qorts/{id}/{id}.QC.summary.txt',
# 	params:
# 		outfile_prefix = '{id}',
# 		outfile_path = '../qorts/{id}/'
# 	log: "../logs/{id}_qorts.log",
# 	benchmark: "../benchmarks/{id}_qorts.benchmark"
# 	shell:
# 		'''
#         java -Xmx20G -jar /home/dvargas/bin/QoRTs.jar QC --singleEnded --outfilePrefix {params.outfile_prefix}. \
#         --skipFunctions writeDESeq --maxReadLength 602 \
#         {input.bam} {input.gtf} {params.outfile_path}  2> {log}
#         '''