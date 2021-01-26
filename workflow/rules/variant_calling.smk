rule lofreq:
    input:
        cram = 'results/data_retrieval/alignment/{accession}.cram',
        crai = 'results/data_retrieval/alignment/{accession}.cram.crai',
        reference = config['reference'],
    output:
        fname_vcf = 'results/variant_calling/calls/{accession}.vcf',
        fname_bam = temp('results/variant_calling/calls/{accession}.bam'),
        fname_bai = temp('results/variant_calling/calls/{accession}.bam.bai'),
    conda:
        '../envs/variant_calling.yaml'
    threads: 8
    shell:
        """
        samtools view \
            -b \
            -T {input.reference} \
            {input.cram} \
        | lofreq viterbi \
            -f {input.reference} \
            - \
        | samtools sort \
            - \
        | lofreq indelqual \
            -f {input.reference} \
            --dindel \
            - \
        > {output.fname_bam}

        # TODO: can I not pipe into this?
        lofreq alnqual \
            -b \
            {output.fname_bam} \
            {input.reference} \
        | sponge {output.fname_bam}

        samtools index {output.fname_bam}

        lofreq call-parallel \
            --pp-threads {threads} \
            -f {input.reference} \
            --call-indels \
            -o {output.fname_vcf} \
            {output.fname_bam}
        """


rule gather_vcf_files:
    input:
        fname_list = expand(
            'results/variant_calling/calls/{accession}.vcf',
            accession=accession_list),
    output:
        fname = 'results/variant_calling/vcf_data.csv'
    script:
        '../scripts/get_base_freqs.py'
