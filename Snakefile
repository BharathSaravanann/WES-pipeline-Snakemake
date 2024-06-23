rule all:
	input:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.analysis.readyGT.funcotated.vcf"


rule bwa_map:
    input:
        ref_file="/home/bharath/Lifecell/ref/hg38.fa",
        fastq1="/home/bharath/Lifecell/fastq/father_R1.fq.gz",
        fastq2="/home/bharath/Lifecell/fastq/father_R2.fq.gz"
    output:
        "output/bwa.sam"
    params:
        sample_id="father",
        rg_id="father",
        platform="ILLUMINA"
    shell:
        """
        bwa mem -t 4 -R "@RG\\tID:{params.rg_id}\\tPL:{params.platform}\\tSM:{params.sample_id}" {input.ref_file} {input.fastq1} {input.fastq2} > {output}
        """
	   
	   
rule sortbam:
	input:
	    "output/bwa.sam"
	    
	    
	output:
	    "output/file.sorted.bam"
	
	
	
	shell:
	    "gatk SortSam -I {input} -O {output} --SORT_ORDER coordinate"
	    
	    
rule indexbam:
	input:
	    "output/file.sorted.bam"
	
	
	
	output:
	    "output/file.sorted.idx"
	
	
	shell:
	    "gatk BuildBamIndex -I {input}"
	    
rule markduplicates:
	input:
	    "output/file.sorted.bam"
	
	
	output:
	    "output/file.sorted.dedup.bam"
	
	
	
	shell:
	    "gatk MarkDuplicatesSpark \
            -I {input} \
            -O {output}"
            
rule BQSR:
	input:
	    bamfile="output/file.sorted.dedup.bam",
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    knownsites="/home/bharath/Lifecell/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
	
	
	output:
	    "output/file.sorted.dedup.recal.data.table"
	
	
	
	shell:
	    "gatk BaseRecalibrator -I {input.bamfile} -R {input.ref_file} --known-sites {input.knownsites} -O {output}"
	  
	  
	  
rule applyBQSR:
	input:
	    bamfile="output/file.sorted.dedup.bam",
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    recal_data="output/file.sorted.dedup.recal.data.table"
	    
	
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.bam"
	
	
	shell:
	    "gatk ApplyBQSR -I {input.bamfile} -R {input.ref_file} --bqsr-recal-file {input.recal_data} -O {output}"
	    
	    
rule Haplotypecaller:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    bamfile="output/file.sorted.dedup.recal.bqsr.bam"
	    
	    
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.vcf"
	
	
	shell:
	    "gatk HaplotypeCaller -R {input.ref_file} -I {input.bamfile} -O {output}"
	    
	   
	   
rule snps:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    vcf_file="output/file.sorted.dedup.recal.bqsr.haplotype.vcf"
	    
	    
	    
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snp.vcf"
	
	
	
	shell:
	    "gatk SelectVariants -R {input.ref_file} -V {input.vcf_file} --select-type SNP -O {output}"
	    
	    
rule indels:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    vcf_file="output/file.sorted.dedup.recal.bqsr.haplotype.vcf"
	    
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.vcf"
	
	
	shell:
	    "gatk SelecVariants -R {input.ref_file} -V {input.vcf_file} --select-type INDEL -O {output}"
	   
	    
	    
	    
rule filtersnps:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    vcf_file="output/file.sorted.dedup.recal.bqsr.haplotype.snp.vcf"
	    
	
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snps.filter.vcf"
	
	
	shell:
	    "gatk VariantFiltration -R {input.ref_file} -V {input.vcf_file} -O {output} \
        -filter-name 'QD_filter' -filter 'QD < 2.0' \
        -filter-name 'FS_filter' -filter 'FS > 60.0' \
        -filter-name 'MQ_filter' -filter 'MQ < 40.0' \
        -filter-name 'SOR_filter' -filter 'SOR > 4.0' \
        -filter-name 'MQRankSum_filter' -filter 'MQRankSum < -12.5' \
        -filter-name 'ReadPosRankSum_filter' -filter 'ReadPosRankSum < -8.0' \
        -genotype-filter-expression 'DP < 10' \
        -genotype-filter-name 'DP_filter' \
        -genotype-filter-expression 'GQ < 10' \
        -genotype-filter-name 'GQ_filter' "
        
     
     
rule filterindels:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    vcf_file="output/file.sorted.dedup.recal.bqsr.haplotype.snp.vcf"
	    
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.vcf"
	
	
	shell:
	    "gatk VariantFiltration -R {input.ref_file} -V {input.vcf_file} -O {output} \
        -filter-name 'QD_filter' -filter 'QD < 2.0' \
        -filter-name 'FS_filter' -filter 'FS > 60.0' \
        -filter-name 'SOR_filter' -filter 'SOR > 4.0' \
        -filter-name 'MQRankSum_filter' -filter 'MQRankSum < -12.5' \
        -genotype-filter-expression 'DP < 10' \
        -genotype-filter-name 'DP_filter' \
        -genotype-filter-expression 'GQ < 10' \
        -genotype-filter-name 'GQ_filter' "
        
	    
rule passsnps:
	input:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snps.filter.vcf"
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snps.filtered.analysis.ready.vcf"
	
	
	
	shell:
	    "gatk SelectVariants \
	--exclude-filtered \
	-V {input} \
	-O {output}"
	
	

rule passindels:
	input:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.vcf"
	    
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.analysis.ready.vcf"
	
	shell:
	    "gatk SelectVariants \
	--exclude-filtered \
	-V {input} \
	-O {output}"
	
	
rule failedgenotypesnps:
	input:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snps.filtered.analysis.ready.vcf"
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snps.filtered.analysis.readyGT.vcf"
	
	
	shell:
	    "cat {input}|grep -v -E 'DP_filter|GQ_filter' > {output}"
	    
	    
	    
rule failedgenotypeindels:
	input:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.analysis.ready.vcf"
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.analysis.readyGT.vcf"
	
	
	shell:
	    "cat {input}|grep -v -E 'DP_filter|GQ_filter' > {output}"
	    
	    
rule snpsannotation:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    vcf_file="output/file.sorted.dedup.recal.bqsr.haplotype.snps.filtered.analysis.readyGT.vcf",
	    data_sources="/home/bharath/germline/funcotator_dataSources.v1.8.hg38.20230908g"
	
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.snps.filtered.analysis.readyGT.funcotated.vcf"
	
	
	
	shell:
	    "gatk Funcotator \
	--variant {input.vcf_file} \
	--reference {input.ref_file} \
	--ref-version hg38 \
	--data-sources-path {input.data_sources} \
	--output {output} \
	--output-file-format VCF"
	    
	    
	    
rule indelsannotation:
	input:
	    ref_file="/home/bharath/Lifecell/ref/hg38.fa",
	    vcf_file="output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.analysis.readyGT.vcf",
	    data_sources="/home/bharath/germline/funcotator_dataSources.v1.8.hg38.20230908g"
	    
	
	
	
	output:
	    "output/file.sorted.dedup.recal.bqsr.haplotype.indels.filtered.analysis.readyGT.funcotated.vcf"
	
	
	shell:
	    "gatk Funcotator \
	--variant {input.vcf_file} \
	--reference {input.ref_file} \
	--ref-version hg38 \
	--data-sources-path {input.data_sources} \
	--output {output} \
	--output-file-format VCF"
	    

	    
	    
	

	    

