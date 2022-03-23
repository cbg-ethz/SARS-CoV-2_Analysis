rule lofreq:
    input:
        cram="results/data_retrieval/alignment/{accession}.cram",
        crai="results/data_retrieval/alignment/{accession}.cram.crai",
        reference=config["reference"],
    output:
        fname_vcf="results/variant_calling/calls/{accession}.vcf",
        fname_bam=temp("results/variant_calling/calls/{accession}.bam"),
        fname_bai=temp("results/variant_calling/calls/{accession}.bam.bai"),
        fname_temp=temp("results/variant_calling/calls/temp/{accession}.another.bam"),
        temp_dir=temp(directory("results/variant_calling/calls/temp_dirs/{accession}/")),
    benchmark:
        "benchmarks/lofreq.{accession}.benchmark.txt"
    conda:
        "../envs/variant_calling.yaml"
    threads: 8
    resources:
        mem_mb=3000,
    shell:
        """
        mkdir -p {output.temp_dir}

        samtools view \
            -b \
            -T {input.reference} \
            --threads {threads} \
            {input.cram} \
        | lofreq viterbi \
            -f {input.reference} \
            - \
        | samtools sort \
            -T {output.temp_dir} \
            --threads {threads} \
            - \
        | lofreq indelqual \
            -f {input.reference} \
            --dindel \
            - \
        > {output.fname_temp}

        # TODO: can I not pipe into this?
        lofreq alnqual \
            -b \
            {output.fname_temp} \
            {input.reference} \
        > {output.fname_bam}

        samtools index -@ {threads} {output.fname_bam}

        lofreq call-parallel \
            --pp-threads {threads} \
            -f {input.reference} \
            --call-indels \
            -o {output.fname_vcf} \
            {output.fname_bam}
        """


rule gather_vcf_files:
    input:
        fname_list=expand(
            "results/variant_calling/calls/{accession}.vcf", accession=accession_list
        ),
        fname_samples="results/data_retrieval/results/selected_samples.csv",
    output:
        fname="results/variant_calling/vcf_data.csv",
    script:
        "../scripts/get_base_freqs.py"
