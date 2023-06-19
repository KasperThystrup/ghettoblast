rule query_copy:
	input:
		fasta	= query_file
	output:
		fasta	= temp("%s/query_idx.fasta" %outdir)
	shell:
		"cp {input.fasta} {output.fasta}"

rule samtools_queryindex:
	input:
		fasta	= rules.query_copy.output.fasta
	output:
		faidx	= temp("%s/query_idx.fasta.fai" %outdir)
	conda:
		"../envs/samtools.yaml"
	message:
		"CMD: samtools faidx {input.fasta}"
	shell:
		"samtools faidx {input.fasta}"
		
