WILDCARDS = glob_wildcards('/ocean/projects/###PATH TO YOUR SAMPLES###/{sample}_R1_001.fastq.gz') #Get list of all sample ids
SALMON_QUANT = expand('03.salmon/{sample}_quant.sf',sample=WILDCARDS.sample) #This will be the final output of the pipeline

rule all: #Tell snakemake to run until salmon outputs are generated 
	input: SALMON_QUANT

rule trim_galore: #Remove adaptors, quality trim, generate fastQC trimming report
	input: 
		r1="/ocean/projects/###PATH TO YOUR SAMPLES###/{sample}_R1_001.fastq.gz",
		r2="/ocean/projects/###PATH TO YOUR SAMPLES###/{sample}_R2_001.fastq.gz"
	output:
		trimmed1="00.qctrim/{sample}_R1_001_val_1.fq.gz",
		trimmed2="00.qctrim/{sample}_R2_001_val_2.fq.gz"
	conda:
		"envs/qc_trim.yaml"
	resources:
		cpus_per_task=4,
		mem_mb=1000
	shell:
		"""
		mkdir -p 00.qctrim/ #make a folder for the outputs to be stored in
		trim_galore -q 20 --cores {resources.cpus_per_task} --phred33 --illumina --length 20 -stringency 3 --fastqc -o 00.qctrim/ --paired {input.r1} {input.r2}
		"""

rule bbduk_rrna: #Remove rRNA from reads
	input:
		trimmed1="00.qctrim/{sample}_R1_001_val_1.fq.gz",
		trimmed2="00.qctrim/{sample}_R2_001_val_2.fq.gz",
		rrna_ref="ref_seqs/smr_v4.3_default_db.fasta"
	output:
		norrna1="01.rrna_remove/{sample}_R1_001_val_1_no_rrna.fq.gz",
		norrna2="01.rrna_remove/{sample}_R2_001_val_2_no_rrna.fq.gz"
	conda:
		"envs/bbtools.yaml"
	resources:
		cpus_per_task=80,
		slurm_partition="RM",
		mem_mb=60000
	shell:
		"""
		mkdir -p 01.rrna_remove/
		bbduk.sh threads={resources.cpus_per_task} in={input.trimmed1} in2={input.trimmed2} k=31 ref={input.rrna_ref} out1={output.norrna1} out2={output.norrna2} minlength=60
		"""

rule bbduk_ercc: #Remove ERCC standards
	input:
		norrna1="01.rrna_remove/{sample}_R1_001_val_1_no_rrna.fq.gz",
		norrna2="01.rrna_remove/{sample}_R2_001_val_2_no_rrna.fq.gz",
		ercc_ref="ref_seqs/ERCC92.fa"
	output:
		noercc1="02.ercc_remove/{sample}_R1_001_val_1_no_rrna_no_ercc.fq.gz",
		noercc2="02.ercc_remove/{sample}_R2_001_val_2_no_rrna_no_ercc.fq.gz"
	conda:
		"envs/bbtools.yaml"
	resources:
		cpus_per_task=20,
		slurm_partition="RM",
		mem_mb=60000
	shell:
		"""
		mkdir -p 02.ercc_remove/
		bbduk.sh threads={resources.cpus_per_task} in={input.norrna1} in2={input.norrna2} k=31 ref={input.ercc_ref} out1={output.noercc1} out2={output.noercc2} minlength=60
		"""

rule salmon_index:
	input: 
		reference_assembly="###PATH TO ASSEMBLY YOU WANT TO MAP TO###"
	output:
		index_bin="03.salmon/assembly_index/pos.bin",
		index_folder=directory("03.salmon/assembly_index")
	conda:
		"envs/salmon.yaml"
	shell:
		"""
		mkdir -p 03.salmon
		salmon index -t {input.reference_assembly} -i {output.index_folder} -k 31 -p 12
		"""

rule salmon_map:
	input:
		index="03.salmon/assembly_index",
		noercc1="02.ercc_remove/{sample}_R1_001_val_1_no_rrna_no_ercc.fq.gz",
		noercc2="02.ercc_remove/{sample}_R2_001_val_2_no_rrna_no_ercc.fq.gz"
	output:
		quant_file = "03.salmon/{sample}_quant.sf"
	conda:
		"salmon"
	resources:
		slurm_partition="RM",
		cpus_per_task=20,
		mem_mb=240000
	shell:
		"""
		salmon quant -i {input.index} -l IU  -1 {input.noercc1} -2 {input.noercc2}  --validateMappings -o 03.salmon/{wildcards.sample}_quant/
		mv 03.salmon/{wildcards.sample}_quant/quant.sf {output.quant_file}
		rm -r 03.salmon/{wildcards.sample}_quant/
		"""
