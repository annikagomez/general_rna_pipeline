WILDCARDS = glob_wildcards('/ocean/projects/ees240006p/agomez3/data/AE1812/{sample}_R1_001.fastq.gz') #Get list of all sample ids
WILDCARDS2 = glob_wildcards('04.merge_drep/merge_drep_95_assembly.{subset}.fasta')
DATABASES = glob_wildcards('databases/{database}.dmnd') #get list of all databases
DIAMOND_ANNOTATIONS = expand('07.annotate/{database}_dmnd.out',database=DATABASES.database)
PHYLODB_ANNOTATIONS = expand("07.annotate/phylodb-estimated-taxonomy_{subset}.out",subset=WILDCARDS2.subset)
EGGNOG_ANNOTATIONS = expand("07.annotate/emapper.annotations.{subset}",subset=WILDCARDS2.subset)
SALMON_QUANT = expand('06.salmon/{sample}_quant.sf',sample=WILDCARDS.sample)

rule all:
	input: PHYLODB_ANNOTATIONS,EGGNOG_ANNOTATIONS,SALMON_QUANT,DIAMOND_ANNOTATIONS

rule trim_galore:
	input: 
		r1="/ocean/projects/ees240006p/agomez3/data/AE1812/{sample}_R1_001.fastq.gz",
		r2="/ocean/projects/ees240006p/agomez3/data/AE1812/{sample}_R2_001.fastq.gz"
	output:
		trimmed1="00.qctrim/{sample}_R1_001_val_1.fq.gz",
		trimmed2="00.qctrim/{sample}_R2_001_val_2.fq.gz"
	conda:
		"qc_trim"
	resources:
		cpus_per_task=4,
		mem_mb=1000
	shell:
		"""
		mkdir -p 00.qctrim/ #make a folder for the outputs to be stored in
		trim_galore -q 20 --cores {resources.cpus_per_task} --phred33 --illumina --length 20 -stringency 3 --fastqc -o 00.qctrim/ --paired {input.r1} {input.r2}
		"""

rule bbduk_rrna:
	input:
		trimmed1="00.qctrim/{sample}_R1_001_val_1.fq.gz",
		trimmed2="00.qctrim/{sample}_R2_001_val_2.fq.gz",
		rrna_ref="ref_seqs/smr_v4.3_default_db.fasta"
	output:
		norrna1="01.rrna_remove/{sample}_R1_001_val_1_no_rrna.fq.gz",
		norrna2="01.rrna_remove/{sample}_R2_001_val_2_no_rrna.fq.gz"
	conda:
		"bbtools"
	resources:
		cpus_per_task=80,
		slurm_partition="RM",
		mem_mb=60000
	shell:
		"""
		mkdir -p 01.rrna_remove/
		bbduk.sh threads={resources.cpus_per_task} in={input.trimmed1} in2={input.trimmed2} k=31 ref={input.rrna_ref} out1={output.norrna1} out2={output.norrna2} minlength=60
		"""

rule bbduk_ercc: 
	input:
		norrna1="01.rrna_remove/{sample}_R1_001_val_1_no_rrna.fq.gz",
		norrna2="01.rrna_remove/{sample}_R2_001_val_2_no_rrna.fq.gz",
		ercc_ref="ref_seqs/ERCC92.fa"
	output:
		noercc1="02.ercc_remove/{sample}_R1_001_val_1_no_rrna_no_ercc.fq.gz",
		noercc2="02.ercc_remove/{sample}_R2_001_val_2_no_rrna_no_ercc.fq.gz"
	conda:
		"bbtools"
	resources:
		cpus_per_task=20,
		slurm_partition="RM",
		mem_mb=60000
	shell:
		"""
		mkdir -p 02.ercc_remove/
		bbduk.sh threads={resources.cpus_per_task} in={input.norrna1} in2={input.norrna2} k=31 ref={input.ercc_ref} out1={output.noercc1} out2={output.noercc2} minlength=60
		"""

rule rnaspades:
	input:
		noercc1="02.ercc_remove/{sample}_R1_001_val_1_no_rrna_no_ercc.fq.gz",
		noercc2="02.ercc_remove/{sample}_R2_001_val_2_no_rrna_no_ercc.fq.gz"
	output:
		assembly="03.rnaspades/{sample}_transcripts.fasta"
	conda:
		"spades"
	resources:
		slurm_partition="RM",
		nodes=1,
		cpus_per_task=18,
		mem_mb=220000
	shell:
		"""
		mkdir -p 03.rnaspades/
		export OMP_NUM_THREADS={resources.cpus_per_task}
		rnaspades.py -1 {input.noercc1} -2 {input.noercc2} -t {resources.cpus_per_task} -o 03.rnaspades/{wildcards.sample}/
		mv 03.rnaspades/{wildcards.sample}/transcripts.fasta {output.assembly}
		rm -r 03.rnaspades/{wildcards.sample}
		"""

rule mmseqs_98:
	input: 
		assembly="03.rnaspades/{sample}_transcripts.fasta"
	output:
		drep_assembly="04.merge_drep/{sample}_drep_98.fasta"
	conda:
		"mmseqs"
	resources:
		slurm_partition="RM-shared",
		nodes=1,
		cpus_per_task=40
	shell:
		"""
		mkdir -p 04.merge_drep
		mmseqs easy-linclust {input.assembly} 04.merge_drep/{wildcards.sample} 04.merge_drep/mmseqs_tmp_{wildcards.sample} --min-seq-id 0.98 --cov-mode 1 -c 0.8 --remove-tmp-files
		mv 04.merge_drep/{wildcards.sample}_rep_seq.fasta {output.drep_assembly}
		rm 04.merge_drep/{wildcards.sample}_all_seqs.fasta
		rm 04.merge_drep/{wildcards.sample}_cluster.tsv
		"""

rule merge_assembly:
	input:
		assemblies=expand("04.merge_drep/{sample}_drep_98.fasta",sample=WILDCARDS.sample)
	output:
		merged_assembly="04.merge_drep/merged_assembly.fasta"
	resources:
		slurm_partition="RM",
		nodes=1,
		cpus_per_task=1,
		mem_mb=200000
	shell:
		"""
		cat {input.assemblies} >> {output.merged_assembly}
		"""

rule mmseqs_db:
	input:
		merged_assembly="04.merge_drep/merged_assembly.fasta"
	output:
		assembly_db="04.merge_drep/merged_assembly_db"
	conda:
		"mmseqs"
	resources:
		slurm_partition="RM",
		nodes=1,
		cpus_per_task=2,
		mem_mb=240000
	shell:
		"mmseqs createdb {input} {output}"

rule mmseqs_cluster:
	input: 
		assembly_db="04.merge_drep/merged_assembly_db"
	params:
		merge_prefix="04.merge_drep/clust_95_db"
	output: 
		assembly_clusters_0="04.merge_drep/clust_95_db.0",
		assembly_clusters_1="04.merge_drep/clust_95_db.1"
	conda:
		"mmseqs"
	resources:
		slurm_partition="RM-512",
		cpus_per_task=2,
		mem_mb=500000
	shell:
		"mmseqs linclust {input} {params.merge_prefix} mmseqs_tmp --min-seq-id 0.95 --cov-mode 1 -c 0.8 --threads 2"

rule mmseqs_write_out:
	input:
		clusters_0="04.merge_drep/clust_95_db.0",
		clusters_1="04.merge_drep/clust_95_db.1",
		assembly_db="04.merge_drep/merged_assembly_db"
	params:
		clusters_prefix="04.merge_drep/clust_95_db"
	output:
		merged_drep_assembly="04.merge_drep/merge_drep_95_assembly.fasta"
	conda:
		"mmseqs"
	resources:
		slurm_partition="RM",
		nodes=1,
		cpus_per_task=10,
		mem_mb=150000
	shell:
		"""
		mmseqs createsubdb {params.clusters_prefix} {input.assembly_db} 04.merge_drep/DB_clu_rep
		mmseqs convert2fasta 04.merge_drep/DB_clu_rep {output}
		"""

rule translate:
	input:
		merged_drep_assembly="04.merge_drep/merge_drep_95_assembly.{subset}.fasta"
	output:
		six_frame_translation="05.translate/{subset}/merged_assembly_6_frame_{subset}.pep"
	conda:
		"emboss"		
	resources:
		slurm_partition="RM",
		cpus_per_task=20,
		mem_mb=240000
	shell:
		"""
		mkdir -p 05.translate
		mkdir -p {wildcards.subset}
		export OMP_NUM_THREADS={resources.cpus_per_task}
		transeq {input.merged_drep_assembly} {output.six_frame_translation} -frame=6
		"""

rule salmon_index:
	input: 
		merged_drep_assembly="04.merge_drep/merge_drep_95_assembly.fasta"
	output:
		index_bin="06.salmon/assembly_index/pos.bin",
		index_folder=directory("06.salmon/assembly_index")
	conda:
		"salmon"
	resources:
                slurm_partition="RM",
                cpus_per_task=12,
                mem_mb=240000
	shell:
		"""
		mkdir -p 06.salmon
		salmon index -t {input.merged_drep_assembly} -i {output.index_folder} -k 31 -p 12
		"""

rule salmon_assembly:
	input:
		index="06.salmon/assembly_index",
		norrna1="01.rrna_remove/{sample}_R1_001_val_1_no_rrna.fq.gz",
		norrna2="01.rrna_remove/{sample}_R2_001_val_2_no_rrna.fq.gz"
	output:
		quant_file = "06.salmon/{sample}_quant.sf"
	conda:
		"salmon"
	resources:
		slurm_partition="RM-512",
		cpus_per_task=20,
		mem_mb=450000
	shell:
		"""
		salmon quant -i {input.index} -l IU  -1 {input.norrna1} -2 {input.norrna2}  --validateMappings -o 06.salmon/{wildcards.sample}_quant/
		mv 06.salmon/{wildcards.sample}_quant/quant.sf {output.quant_file}
		rm -r 06.salmon/{wildcards.sample}_quant/
		"""

rule eukulele:
	input:
		six_frame_translation="05.translate/{subset}/merged_assembly_6_frame_{subset}.pep"
	output:
		tax_est="07.annotate/phylodb-estimated-taxonomy_{subset}.out"
	conda:
		"new_EUKulele"
	resources:
		slurm_partition="RM-shared",
		nodes=1,
		cpus_per_task=30,
		mem_mb=48000
	shell:
		"""
		mkdir -p 07.annotate
		EUKulele --sample_dir 05.translate/{wildcards.subset} -m mets --protein_extension .pep --database phylodb --reference_dir databases/phylodb -o 07.annotate/eukulele_out_{wildcards.subset} --no_busco
		mv 07.annotate/eukulele_out_{wildcards.subset}/taxonomy_estimation/merged_assembly_6_frame_{wildcards.subset}-estimated-taxonomy.out {output.tax_est}
		"""

rule eggnog:
	input: 
		six_frame_translation="05.translate/{subset}/merged_assembly_6_frame_{subset}.pep"
	output:
		annotations="07.annotate/emapper.annotations.{subset}"
	conda:
		"eggnog"
	resources:
		slurm_partition="RM-shared",
		nodes=1,
		cpus_per_task=30,
		mem_mb=48000
	shell:
		"""
		mkdir -p 07.annotate/eggnog_out_{wildcards.subset}
		emapper.py -i {input.six_frame_translation} --itype proteins -o split --output_dir 07.annotate/eggnog_out_{wildcards.subset}/ --cpu 30
		mv 07.annotate/eggnog_out_{wildcards.subset}/split.emapper.annotations {output.annotations}
		"""

rule diamond_blastp:
	input:
		six_frame_translation="05.translate/merged_assembly_6_frame.pep",
		database="databases/{database}.dmnd"		
	output:
		diamond_results="07.annotate/{database}_dmnd.out"
	conda:
		"diamond"
	resources:
		slurm_partition="RM",
		nodes=1,
		cpus_per_task=64
	shell:
		"""
		mkdir -p 07.annotate
		diamond blastp -d {input.database} -q {input.six_frame_translation} -o {output.diamond_results}
		"""
