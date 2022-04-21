"""
Query sequences are blasted against each individual sample fasta file.
Results are saved to a custom build results table for all samplese.
"""
rule blast_all_fasta:
	input:
		fasta	= sample_path+"{sample}.fasta"
	params:
		query	= config["query_file"]
	output:
		blast	= temp(outdir+"{sample}.blast")
	conda:
		"../envs/blast.yaml"
	shell:
		"blastn -query {params.query} -subject {input.fasta} -out {output.blast} -outfmt '6 pident nident length mismatch gapopen bitscore qacc qstart qend qseq sacc sstart send sseq'"


rule blast_all_fna:
	input:
		fasta	= sample_path+"{sample}.fna"
	params:
		query	= config["query_file"]
	output:
		blast	= temp(outdir+"{sample}.blast")
	conda:
		"../envs/blast.yaml"
	shell:
		"blastn -query {params.query} -subject {input.fasta} -out {output.blast} -outfmt '6 pident nident length mismatch gapopen bitscore qacc qstart qend qseq sacc sstart send sseq'"


"""
Rename each blast.
"""
rule save_blast:
	input:
		blast = outdir+"{sample}.blast"
	output:
		blast_out = outdir+"{sample}_blast.tsv"
	shell:
		"mv {input.blast} {output.blast_out}"

"""
Only perform blast and save generated blast files.
"""
rule blast:
	input:
		expand(rules.save_blast.output.blast_out, sample = samples)

"""
Delete saved blast files
"""
rule delete_saved_blast:
	input:
		delete = rules.save_blast.output.blast_out
	shell:
		"rm {input}"

"""
Clear up all saved blast files
"""
rule clear_blast:
	input:
		expand(rules.delete_saved_blast.input.delete, sample = samples)
