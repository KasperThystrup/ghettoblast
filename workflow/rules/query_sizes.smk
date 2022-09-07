rule samtools_queryindex:
	input:
		fasta	= query_file
	output:
		faidx	= temp("%s/query_idx.fai" %outdir)
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools faidx {input.fasta} --output {output.faidx}"
		