rule_name = "download_sequences_via_lapis"
rule download_sequences_via_lapis:
    """Download sequences via LAPIS"""

    message: """Downloading sequences via LAPIS\n

    Log:           {log}

    Output:
        sequences: {output.sequences}
    """   

    output:
        sequences = "data/{build_name}_sequences.fasta",
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        curl https://mpox-lapis.genspectrum.org/v1/sample/fasta --output {output.sequences}
        """

rule_name = "download_metadata_via_lapis"
rule download_metadata_via_lapis:
    """Download metadata via LAPIS"""

    message: """Downloading metadata via LAPIS\n

    Log:          {log}

    Output:
        metadata: {output.metadata}
    """ 

    output:
        metadata = "data/{build_name}_metadata.tsv"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        curl https://mpox-lapis.genspectrum.org/v1/sample/details?dataFormat=csv | \
            tr -d "\r" |
            sed -E 's/("([^"]*)")?,/\\2\\t/g' > {output.metadata}
        """
