"""
Query sequences are blasted against each individual sample fasta file.
Results are saved to a custom build results table for all samplese.
"""
rule blast_all:
	input:
		fasta	= "%s/{sample}.fasta" %sample_path
	params:
		query	= config["query_file"]
	output:
		blast	= temp("%s/tmp/{sample}.blast" %outdir)
	conda:
		"../envs/blast.yaml"
	message:
		"CMD: blastn -query {params.query} -subject {input.fasta} -out {output.blast} -outfmt '6 pident nident length mismatch gapopen bitscore qacc qstart qend qseq sacc sstart send sseq' %s" %besthit
	shell:
		"blastn -query {params.query} -subject {input.fasta} -out {output.blast} -outfmt '6 pident nident length mismatch gapopen bitscore qacc qstart qend qseq sacc sstart send sseq' %s" %besthit
