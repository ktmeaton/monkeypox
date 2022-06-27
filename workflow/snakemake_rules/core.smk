'''
This part of the workflow expects input files

        sequences = "data/sequences.fasta",
        metadata =  "data/metadata.tsv",

and will produce output files as

        auspice_json = auspice_dir + "/monkeypox_{build_name}.json"

Parameter are expected to sit in the `config` data structure.
In addition, `build_dir` and `auspice_dir` need to be defined upstream.
'''

# -----------------------------------------------------------------------------
rule_name = "wrangle_metadata"
rule wrangle_metadata:
    """Wrangle metadata."""
    message: """Wrangling metadata.

    Log:                     {log}

    Input:
        metadata:            {input.metadata}

    Output:
        metadata:            {output.metadata}
    """
    input:
        metadata =  "data/{build_name}_metadata.tsv"
    output:
        metadata = "results/{build_name}_metadata.tsv"
    params:
        strain_id = lambda w: config.get('strain_id_field', 'strain')
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        python3 scripts/wrangle_metadata.py \
            --metadata {input.metadata} \
            --strain-id {params.strain_id} \
            --output {output.metadata} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "filter"
rule filter:
    """Filter metadata and sequences."""
    message: """Filtering metadata and sequences.

    Log:                     {log}

    Input:
        sequences:           {input.sequences}
        metadata:            {input.metadata}
        exclude:             {input.exclude}

    Output:
        sequences:           {output.sequences}
        metadata:            {output.metadata}
        log:                 {output.log}

    Params:
        group_by:            {params.group_by!s}
        sequences_per_group: {params.sequences_per_group}
        min_date:            {params.min_date}
        min_length:          {params.min_length}
    """
    input:
        sequences = "data/{build_name}_sequences.fasta",
        metadata =  "results/{build_name}_metadata.tsv",
        exclude = config["exclude"]
    output:
        sequences = build_dir + "/{build_name}/filtered.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv",
        log = build_dir + "/{build_name}/filtered.log",
    params:
        group_by = "country year",
        sequences_per_group = 1000,
        min_date = config['min_date'],
        min_length = config['min_length']
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --output-log {output.log} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "align"
rule align:
    """Align sequences to reference."""
    message: """Aligning sequences to reference and filling gaps with N.

    Log:              {log}

    Input:
        sequences:    {input.sequences}
        reference:    {input.reference}

    Output:
        alignment:    {output.alignment}
        insertions:   {output.insertions}

    Params:
        max_indel:    {params.max_indel}
        seed_spacing: {params.seed_spacing}
    """
    input:
        sequences = rules.filter.output.sequences,
        reference = config["reference"]
    output:
        alignment = build_dir + "/{build_name}/aligned.fasta",
        insertions = build_dir + "/{build_name}/insertions.fasta"
    params:
        max_indel = config["max_indel"],
        seed_spacing = config["seed_spacing"]
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        nextalign2 run \
            --jobs {resources.cpus} \
            --reference {input.reference} \
            --max-indel {params.max_indel} \
            --seed-spacing {params.seed_spacing} \
            --output-fasta {output.alignment} \
            --output-insertions {output.insertions} \
            {input.sequences} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "mask"
rule mask:
    """Mask ends of the alignment"""
    message: """Masking ends of the alignment.

    Log:              {log}

    Input:
        sequences:    {input.sequences}
        mask:         {input.mask}

    Output:
        sequences:    {output.sequences}

    Params:
        from_start:   {params.from_start}
        from_end:     {params.from_end}
    """
    input:
        sequences = rules.align.output.alignment,
        mask = config["mask"]["maskfile"]
    output:
        sequences = build_dir + "/{build_name}/masked.fasta"
    params:
        from_start = config["mask"]["from_beginning"],
        from_end = config["mask"]["from_end"]
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur mask \
          --sequences {input.sequences} \
          --mask {input.mask} \
          --mask-from-beginning {params.from_start} \
          --mask-from-end {params.from_end} \
          --output {output.sequences} \
          > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "tree"

def _params_tree():
    ""
    params = {
        "tree_builder_args" : "",
        "override_default_args" : ""
    }
    if "tree_builder_args" in config:
        tree_builder_args = config["tree_builder_args"]
        if tree_builder_args:
            params["tree_builder_args"] = "--tree-builder-args=\"{}\"".format(tree_builder_args)

    if "override_default_args" in config:
        override_default_args = config["override_default_args"]
        if override_default_args == True:
            params["override_default_args"] = "--override-default-args"

    return params


rule tree:
    """Build tree"""
    message: """Building tree.
    Log:                       {log}

    Input:
        alignment:             {input.alignment}

    Output:
        tree:                  {output.tree}

    Params:
        override_default_args: {params.override_default_args}    
        tree_builder_args:     {params.tree_builder_args}
    """
    input:
        alignment = rules.mask.output.sequences,
    output:
        tree = build_dir + "/{build_name}/tree_raw.nwk",
    params:
        tree_builder_args = lambda w: _params_tree()["tree_builder_args"],
        override_default_args = lambda w: _params_tree()["override_default_args"],
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {resources.cpus} \
            {params.override_default_args} \
            {params.tree_builder_args} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "refine"
rule refine:
    """Refine tree to estimate a timetree."""
    message: """Refining tree to estimate a timetree.

    Log:                  {log}

    Input:
        alignment:        {input.alignment}
        metadata:         {input.metadata}
        tree:             {input.tree}

    Output:
        tree:             {output.tree}
        node_data:        {output.node_data}

    Params:
        coalescent:       {params.coalescent}
        clock_filter_iqd: {params.clock_filter_iqd}
        clock_rate:       {params.clock_rate}
        clock_std_dev:    {params.clock_std_dev}
        date_inference:   {params.date_inference}
        root:             {params.root}
    """
    input:
        tree = rules.tree.output.tree,
        alignment = build_dir + "/{build_name}/masked.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        tree = build_dir + "/{build_name}/tree.nwk",
        node_data = build_dir + "/{build_name}/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 0,
        root = config["root"],
        clock_rate = lambda w: f"--clock-rate {config['clock_rate']}" if "clock_rate" in config else "",
        clock_std_dev = lambda w: f"--clock-std-dev {config['clock_std_dev']}" if "clock_std_dev" in config else ""
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --timetree \
            --root {params.root} \
            --precision 3 \
            --keep-polytomies \
            {params.clock_rate} \
            {params.clock_std_dev} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --date-confidence \
            --clock-filter-iqd {params.clock_filter_iqd} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "ancestral"
rule ancestral:
    """Reconstruct ancestral sequences and mutations."""
    message: """Reconstructing ancestral sequences and mutations.

    Log:           {log}

    Input:
        alignment: {input.alignment}
        tree:      {input.tree}

    Output:
        node_data: {output.node_data}

    Params:
        inference: {params.inference}
    """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output.sequences,
    output:
        node_data = build_dir + "/{build_name}/nt_muts.json"
    params:
        inference = "joint"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "translate"
rule translate:
    """Translate amino acid sequences."""
    message: """Translating amino acid sequences.

    Log:                  {log}

    Input:
        tree:              {input.tree}
        node_data:         {input.node_data}
        genbank_reference: {input.genbank_reference}

    Output:
        node_data:        {output.node_data}
    """
    input:
        tree              = rules.refine.output.tree,
        node_data         = rules.ancestral.output.node_data,
        genbank_reference = config["genbank_reference"]
    output:
        node_data         = build_dir + "/{build_name}/aa_muts.json"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.genbank_reference} \
            --output {output.node_data} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "traits"
rule traits:
    """Infer ancestral traits."""
    message: """Inferring ancestral traits.

    Log:                          {log}

    Input:
        tree:                     {input.tree}
        metadata:                 {input.metadata}

    Output:
        node_data:                {output.node_data}

    Params:
        columns:                  {params.columns}
        sampling_bias_correction: {params.sampling_bias_correction}
    """
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        node_data = build_dir + "/{build_name}/traits.json",
    params:
        columns = "country",
        sampling_bias_correction = 3
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction}
        """

# -----------------------------------------------------------------------------
rule_name = "clades"
rule clades:
    """Add internal clade labels."""
    message: """Adding internal clade labels.

    Log:           {log}

    Input:
        tree :     {input.tree}
        aa_muts:   {input.aa_muts}
        nuc_muts:  {input.nuc_muts}
        clades:    {input.clades}

    Output:
        node_data: {output.node_data}
    """
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = config["clades"]
    output:
        node_data = build_dir + "/{build_name}/clades.json"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "mutation_context"
rule mutation_context:
    """Generate mutation context"""
    message: """Generating mutation context.

    Log:           {log}

    Input:
        tree:      {input.tree}
        node_data: {input.node_data}

    Output:
        node_data: {output.node_data}
    """
    input:
        tree = rules.refine.output.tree,
        node_data = build_dir + "/{build_name}/nt_muts.json"
    output:
        node_data = build_dir + "/{build_name}/mutation_context.json",
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        python3 scripts/mutation_context.py \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --output {output.node_data} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "remove_time"
rule remove_time:
    """Remove time information"""
    message: """Removing time information.

    Log:           {log}

    Input:
        tree:      {input.tree}
        node_data: {input.node_data}

    Output:
        node_data: {output.node_data}
    """
    input:
        node_data = "results/{build_name}/branch_lengths.json"
    output:
        node_data = "results/{build_name}/branch_lengths_no_time.json"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        python3 scripts/remove_timeinfo.py \
          --input-node-data {input.node_data} \
          --output-node-data {output.node_data} \
          > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "recency"
rule recency:
    """Construct submission recency field."""
    message: """Constructing submission recency field.

    Log:           {log}

    Input:
        metadata:  {input.metadata}

    Output:
        node_data: {output.node_data}
    """
    input:
        metadata = build_dir + "/{build_name}/metadata.tsv"
    output:
        node_data = build_dir + "/{build_name}/recency.json"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output.node_data} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "export"
rule export:
    """Export data files for for auspice"""
    message: """Exporting data files for for auspice.

    Log:                  {log}

    Input:
        aa_muts:          {input.aa_muts}
        auspice_config:   {input.auspice_config}
        branch_lengths:   {input.branch_lengths}
        clades:           {input.clades}
        colors:           {input.colors}
        description:      {input.description}
        lat_longs:        {input.lat_longs}
        metadata:         {input.metadata}
        mutation_context: {input.mutation_context}
        nt_muts:          {input.nt_muts}
        recency:          {input.recency}
        traits:           {input.traits}
        tree:             {input.tree}

    Output:
        auspice_json:     {output.auspice_json}
        root_sequence:    {output.root_sequence}
    """
    input:
        tree = rules.refine.output.tree,
        metadata = build_dir + "/{build_name}/metadata.tsv",
        branch_lengths = lambda w: "results/{build_name}/branch_lengths.json" if config.get('timetree', False) else "results/{build_name}/branch_lengths_no_time.json",
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.node_data,
        mutation_context = rules.mutation_context.output.node_data,
        recency = lambda w: rules.recency.output.node_data if config.get('recency', False) else [],
        colors = config["colors"],
        lat_longs = config["lat_longs"],
        description = config["description"],
        auspice_config = config["auspice_config"]
    output:
        auspice_json =  build_dir + "/{build_name}/raw_tree.json",
        root_sequence = build_dir + "/{build_name}/raw_tree_root-sequence.json"
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.mutation_context} {input.clades} {input.recency} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json} \
            > {log} 2>&1;
        """

# -----------------------------------------------------------------------------
rule_name = "final_strain_name"
rule final_strain_name:
    """Construct final strain name."""
    message: """Constructing submission recency field.
    Log:           {log}

    Input:
        auspice_json:     {input.auspice_json}
        metadata:         {input.metadata}
        root_sequence:    {input.root_sequence}

    Output:
        auspice_json:     {output.auspice_json}
        root_sequence:    {output.root_sequence}
    """

    input:
        auspice_json =  build_dir + "/{build_name}/raw_tree.json",
        metadata = build_dir + "/{build_name}/metadata.tsv",
        root_sequence = build_dir + "/{build_name}/raw_tree_root-sequence.json"
    output:
        auspice_json =  build_dir + "/{build_name}/tree.json",
        root_sequence =  build_dir + "/{build_name}/tree_root-sequence.json"
    params:
        display_strain_field = lambda w: config.get('display_strain_field', 'strain')
    benchmark:
        "benchmarks/{rule}/{{build_name}}_{today}.tsv".format(today=today, rule=rule_name),
    log:
        "logs/{rule}/{{build_name}}_{today}.log".format(today=today, rule=rule_name),
    shell:
        """
        python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
                --input-auspice-json {input.auspice_json} \
                --display-strain-name {params.display_strain_field} \
                --output {output.auspice_json} \
                > {log} 2>&1;
        cp {input.root_sequence} {output.root_sequence}
        """
