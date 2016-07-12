import os

# Dependencies

dependencies = config['dependencies']

trim_galore= dependencies['trim_galore']
cutadapt = dependencies['cutadapt']
sambamba = dependencies['sambamba']
markDup = dependencies['markDup']
gatk = dependencies['gatk']
genomeFile = dependencies['genomeFile']a
star = dependencies['star']
star_genome_dir  = dependencies['star_genome_dir']
## References

references = config['references']

referenceIndel = references['referenceIndel']
referenceSnp = references['referenceSnp']

samples = config['samples']
cwd = os.getcwd()

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
   






