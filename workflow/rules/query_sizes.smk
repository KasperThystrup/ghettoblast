rule samtools_queryindex:
	input:
		fasta	= query_file
	output:
		faidx	= temp(outdir+"query_idx.fai")
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools faidx {input.fasta} --output {output.faidx}"

rule queryindex:
	input:
		query_index	= outdir+"query_idx.fai"
