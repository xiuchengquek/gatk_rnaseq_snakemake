import os

# Dependencies
trim_galore = "/share/ClusterShare/software/contrib/gi/trim_galore/0.4.0/trim_galore"
cutadapt = "/share/ClusterShare/software/contrib/fabbus/python/2.7.3/bin/cutadapt"
sambamba = "/home/xiuque/bin/sambamba_v0.5.8"
markDup = "/home/xiuque/bin/MarkDuplicates.jar"
gatk = "/home/xiuque/bin/GenomeAnalysisTK.jar"
genomeFile = "/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.8/b37/human_g1k_v37.fasta"
referenceIndel = "/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
referenceSnp = "/share/ClusterShare/biodata/contrib/gi/gatk-resource-bundle/2.8/b37/dbsnp_138.b37.vcf"

cwd = os.getcwd()

star = "/share/ClusterScratch/xiuque/testsnakemake/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR"
star_genome_dir  = "/share/ClusterScratch/xiuque/testsnakemake/lib/genomedir/"

samples = ["778-parental-04-02-2014","778-parental-08-04-2014","778-parental-16-05-2014","778-parental-24-03-2014","778-R-CDK4i-03-01-2014","778-R-CDK4i-04-02-2014","778-R-CDK4i-08-07-2014","778-R-CDK4i-16-05-2014","778-R-CDK4i-24-03-2014","778-R-Nutlin-03-01-2014","778-R-Nutlin-04-02-2014","778-R-Nutlin-08-07-2014","778-R-Nutlin-16-05-2014","778-R-Nutlin-24-03-2014","778-R-Tunicamycin-03-01-2014","778-R-Tunicamycin-04-02-2014","778-R-Tunicamycin-08-07-2014","778-R-Tunicamycin-16-05-2014","778-R-Tunicamycin-24-03-2014"]

rule targets:
    input  : expand("data/{sample}/trimmed/Aligned.recali.bam", sample=samples) , expand("data/{sample}/trimmed/Aligned.recaliNoDS.bam", sample=samples)
   
rule trimgalore:
    input  : "data/{sample, [\w\-]+}/{sample}_R1.fastq.gz", "data/{sample, [\w\-]+}/{sample}_R2.fastq.gz", "data/{sample, [\w\-]+}/"
    output  : "data/{sample, [\w\-]+}/trimmed/{sample}_R1_val_1.fq.gz", "data/{sample, [\w\-]+}/trimmed/{sample}_R2_val_2.fq.gz"
    params: core="2", mem="8g"
    shell  :  "{trim_galore} --paired --gzip --path_to_cutadapt {cutadapt} --o {input[2]}/trimmed {input[0]} {input[1]}"  

rule star:
    input : "data/{sample, [\w\-]+}/trimmed/{sample}_R1_val_1.fq.gz", "data/{sample, [\w\-]+}/trimmed/{sample}_R2_val_2.fq.gz", "data/{sample, [\w\-]+}"
    output : "data/{sample, [\w\-]+}/trimmed/Aligned.out.bam"
    params: core="16", sampleName = os.path.basename("data/{sample, [\w\-]+}"), wd = "data/{sample, [\w\-]+}/trimmed"
    log : "star.log"
    shell : "cd {params.wd} && {star} --genomeDir {star_genome_dir} --readFilesCommand zcat --readFilesIn {cwd}/{input[0]} {cwd}/{input[1]} --runThreadN {params.core}" +		" --outSAMattrRGline ID:{params.sampleName} SM:{params.sampleName} PL:ILLUMINA PU:HISEQ2000 LB:GI --outSAMtype BAM Unsorted;"

rule sortBam:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.out.bam"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.sorted.bam"   
    params: core="8", sampleName = os.path.basename("data/{sample, [\w\-]+}"), wd = "data/{sample, [\w\-]+}/trimmed"
    shell: "{sambamba} sort --memory-limit=60G --out={output} --nthreads={params.core} {input}"

rule markDup:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.sorted.bam"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.marked.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.mark.bai"
    params: core="4" , wd = "data/{sample, [\w\-]+}/trimmed/"
    shell : "java -jar -Xmx24G {markDup} I={input} O={output[0]} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={params.wd}/output.metrics" 

rule splitNCigarReads:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.marked.bam"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.split.bam" , "data/{sample, [\w\-]+}/trimmed/Aligned.split.bai"
    params: core="2" , wd = "data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar -Xmx8G {gatk} -T SplitNCigarReads -R {genomeFile} -I {input} -o {output[0]} -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60"    

rule makeIntervalFile:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.split.bam"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.target.intervals"
    params: core="8", wd = "data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar -Xmx56G {gatk} -T RealignerTargetCreator -R {genomeFile} -I {input} --known {referenceIndel} -o {output} -nt 16" 

rule indelRealigner:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.split.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.target.intervals"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.indel.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.indel.bai"
    params: core="2", wd = "data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar -Xmx8G {gatk} -T IndelRealigner -R {genomeFile} -I {input[0]} -known {referenceIndel} -targetIntervals {input[1]} -o {output[0]}"

rule baseRecalibrator: 
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.indel.bam"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.recali.table"
    params: core="8", wd = "data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar -Xmx56G {gatk} -T BaseRecalibrator -R {genomeFile} -I {input} -L 20 -knownSites {referenceSnp} -knownSites {referenceIndel} -o {output} -nct 8"

rule analyzeCovariates:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.recali.table"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.recali.post.table", "data/{sample, [\w\-]+}/trimmed/recalibration_plots.pdf"
    params: core="2", wd = "data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar AnalyzeCovariates -R {genomeFile} -L 20 -before {input} -after {output[0]} -plots {output[1]}"

rule printReads:
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.indel.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.recali.table"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.recali.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.recali.bai"
    params: core="8" , wd = "data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar -Xmx56G {gatk} -T PrintReads -R {genomeFile} -I {input[0]} -o {output[0]} -L 20 -BQSR {input[1]} -nct 8"
    
rule printReadsNoDownSampling:   
    input: "data/{sample, [\w\-]+}/trimmed/Aligned.indel.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.recali.table"
    output: "data/{sample, [\w\-]+}/trimmed/Aligned.recaliNoDS.bam", "data/{sample, [\w\-]+}/trimmed/Aligned.recaliNoDS.bai"
    params: core="8", wd="data/{sample, [\w\-]+}/trimmed/"
    shell: "java -jar -Xmx56G {gatk} -T PrintReads -R {genomeFile} -I {input[0]} -o {output[0]} -BQSR {input[1]} -nct 8"
   






