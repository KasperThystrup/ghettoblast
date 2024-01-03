rule merge_blast:
	input:
		blast		= expand(rules.blast_all.output.blast, sample = samples),
		faidx		= rules.samtools_queryindex.output.faidx
	params:
		metadata	= metadata,
		exclude_seqs	= exclude_seqs,
		top_only	= top_only,
		outdir		= outdir,
		debug		= debug
	output:
		blast_results	= "%s/%s_blast.tsv" %(outdir, query_name)
	conda:
		"../envs/R.yaml"
	script:
		"../scripts/merge_blast.R"
