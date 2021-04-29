rule isolate_AD:
    input:
        fa=config['common']['genome'],
    output:
        fa="resources/chrM.fa",
        idx="resources/chrM.fa.fai",
    conda:
        "../envs/samtools.yaml",
    shell:
        "samtools faidx {input.fa} chrM > {output.fa}; "
        "samtools faidx {output.fa}; "
