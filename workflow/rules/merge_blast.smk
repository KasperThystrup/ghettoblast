rule merge_blast:
	input:
		blast		= expand(outdir+"{sample}.blast", sample = samples),
		faidx		= rules.samtools_queryindex.output.faidx
	params:
		metadata	= metadata
	output:
		blast_results	= outdir+query_name+"_blast.tsv"
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/merge_blast.R"
	
rule merge:
	input:
		blast		= rules.blast.input
	params:
		blast_columns   = "config/blast_columns.tsv",
		metadata        = metadata
	output:
		blast_results   = outdir+query_name+"__blast.tsv"
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/merge_blast.R"

