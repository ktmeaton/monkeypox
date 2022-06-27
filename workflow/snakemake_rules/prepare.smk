rule_name = "download"
rule download:
    """Download sequences and metadata from data.nextstrain.org"""

    message: """Downloading sequences and metadata from data.nextstrain.org\n

    Log:               {log}

    Output:
        sequences:     {output.sequences}
        metadata:      {output.metadata}
    
    Params:
        sequences_url: {params.sequences_url}
        metadata_url:  {params.metadata_url}
        strain_id:     {params.strain_id}
    """    
    output:
        sequences = "data/{build_name}_sequences.fasta.xz",
        metadata = "data/{build_name}_metadata.tsv.gz",
    params:
        sequences_url = "https://data.nextstrain.org/files/workflows/monkeypox/sequences.fasta.xz",
        metadata_url = "https://data.nextstrain.org/files/workflows/monkeypox/metadata.tsv.gz",
        strain_id = lambda w: config.get('strain_id_field', 'strain'),
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        curl -fsSL --compressed {params.sequences_url:q} --output {output.sequences} > {log} 2>&1;
        curl -fsSL --compressed {params.metadata_url:q} --output {output.metadata} >> {log} 2>&1;
        """

rule_name = "decompress"
rule decompress:
    """Decompress sequences and metadata"""

    message: """Decompressing sequences and metadata\n

    Log:               {log}

    Input:
        sequences:     {input.sequences}
        metadata:      {input.metadata}

    Output:
        sequences:     {output.sequences}
        metadata:      {output.metadata}
    """    


    input:
        sequences = "data/{build_name}_sequences.fasta.xz",
        metadata = "data/{build_name}_metadata.tsv.gz",
    output:
        sequences = "data/{build_name}_sequences.fasta",
        metadata = "data/{build_name}_metadata.tsv",
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        gzip --decompress --keep {input.metadata} > {log} 2>&1;
        xz -q -Q --decompress --keep {input.sequences} >> {log} 2>&1;
        """
