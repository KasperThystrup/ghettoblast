rule merge_blast:
	input:
		blast		= expand(rules.blast_all.output.blast, sample = samples),
		faidx		= rules.samtools_queryindex.output.faidx
	params:
		blast_columns	= "config/blast_columns.tsv",
		metadata	= metadata
	output:
		blast_results	= outdir+query_name+"_blast.tsv"
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/merge_blast.R"
	
#rule merge:
#	input:
#		rules.save_blast.output.blast_out
#	params:
#                blast_columns   = "config/blast_columns.tsv",
#                metadata        = metadata
#        output:
#                blast_results   = outdir+query_name+"_blast.tsv"
#        conda:
#                "../envs/R.yaml"
#        script:
#                "../scripts/merge_blast.R"

