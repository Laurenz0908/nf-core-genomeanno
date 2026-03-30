process DESIGN_STRAIN_PRIMERS {
    tag "${meta.id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/pip_biopython_primer3-py:e339c879b0de559a'

    input:
    tuple val(meta), path(target), path(bed_files)

    output:
    tuple val(meta), path("${prefix}_primers.tsv")  , emit: primers_tsv
    tuple val(meta), path("${prefix}_primers.json")  , emit: primers_json
    tuple val(meta), path("${prefix}_primers.fasta") , emit: primers_fasta , optional: true
    tuple val(meta), path("${prefix}_unique.bed")    , emit: unique_bed
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def min_unique_length = task.ext.min_unique_length ?: 200
    def primer_opt_tm     = task.ext.primer_opt_tm     ?: 60.0
    def primer_min_tm     = task.ext.primer_min_tm     ?: 57.0
    def primer_max_tm     = task.ext.primer_max_tm     ?: 63.0
    def primer_min_gc     = task.ext.primer_min_gc     ?: 40.0
    def primer_max_gc     = task.ext.primer_max_gc     ?: 60.0
    def product_size_min  = task.ext.product_size_min  ?: 150
    def product_size_max  = task.ext.product_size_max  ?: 500
    def num_primer_pairs  = task.ext.num_primer_pairs  ?: 3
    """
    #!/usr/bin/env python3
    import csv
    import json
    import sys
    from dataclasses import dataclass, asdict
    from pathlib import Path
    from subprocess import run

    import primer3
    from Bio import SeqIO

    prefix             = "${prefix}"
    target_fasta       = "${target}"
    min_unique_length  = ${min_unique_length}

    # ── Step 1: Merge all BED files ──────────────────────────────────────────
    bed_files = sorted(Path(".").glob("*.bed"))
    if not bed_files:
        sys.exit("ERROR: No BED files found")

    # Concatenate, sort, merge
    with open("combined.bed", "w") as out:
        for bf in bed_files:
            out.write(bf.read_text())

    run("sort -k1,1 -k2,2n combined.bed > sorted.bed", shell=True, check=True)
    run("cat sorted.bed | python3 -c '"
        "import sys\\n"
        "prev = None\\n"
        "for line in sys.stdin:\\n"
        "    c,s,e = line.strip().split(chr(9))[:3]\\n"
        "    s,e = int(s),int(e)\\n"
        "    if prev and prev[0]==c and s<=prev[2]:\\n"
        "        prev = (c, prev[1], max(e,prev[2]))\\n"
        "    else:\\n"
        "        if prev: print(f\"{prev[0]}\\\\t{prev[1]}\\\\t{prev[2]}\")\\n"
        "        prev = (c,s,e)\\n"
        "if prev: print(f\"{prev[0]}\\\\t{prev[1]}\\\\t{prev[2]}\")\\n"
        "' > merged.bed", shell=True, check=True)

    # ── Step 2: Build genome file and find complement ────────────────────────
    seqs = {}
    genome_entries = []
    for rec in SeqIO.parse(target_fasta, "fasta"):
        seqs[rec.id] = str(rec.seq)
        genome_entries.append((rec.id, len(rec.seq)))

    # Sort genome file same order as the sequences
    with open("genome.txt", "w") as fh:
        for name, length in sorted(genome_entries):
            fh.write(f"{name}\\t{length}\\n")

    # Sort merged.bed to match genome.txt ordering
    run("sort -k1,1 -k2,2n merged.bed > merged_sorted.bed", shell=True, check=True)

    # Compute complement (unique regions)
    # Pure Python implementation to avoid bedtools dependency in this container
    covered = {}
    with open("merged_sorted.bed") as fh:
        for line in fh:
            parts = line.strip().split("\\t")
            c, s, e = parts[0], int(parts[1]), int(parts[2])
            covered.setdefault(c, []).append((s, e))

    unique_regions = []
    with open(f"{prefix}_unique.bed", "w") as fh:
        for name, length in sorted(genome_entries):
            intervals = sorted(covered.get(name, []))
            pos = 0
            for s, e in intervals:
                if s > pos and (s - pos) >= min_unique_length:
                    fh.write(f"{name}\\t{pos}\\t{s}\\n")
                    unique_regions.append((name, pos, s))
                pos = max(pos, e)
            if length > pos and (length - pos) >= min_unique_length:
                fh.write(f"{name}\\t{pos}\\t{length}\\n")
                unique_regions.append((name, pos, length))

    print(f"Found {len(unique_regions)} unique regions >= {min_unique_length} bp")

    # ── Step 3: Design primers ───────────────────────────────────────────────
    @dataclass
    class PrimerPair:
        region_id: str; contig: str; region_start: int; region_end: int
        fwd_seq: str; rev_seq: str; fwd_tm: float; rev_tm: float
        fwd_gc: float; rev_gc: float; product_size: int; penalty: float

    primers = []
    product_size_max = ${product_size_max}

    for idx, (contig, start, end) in enumerate(unique_regions, 1):
        rid = f"region_{idx:04d}"
        region_seq = seqs[contig][start:end]

        max_template = product_size_max + 200
        if len(region_seq) > max_template:
            offset = (len(region_seq) - max_template) // 2
            region_seq = region_seq[offset:offset + max_template]

        try:
            result = primer3.design_primers(
                seq_args={"SEQUENCE_ID": rid, "SEQUENCE_TEMPLATE": region_seq},
                global_args={
                    "PRIMER_NUM_RETURN": ${num_primer_pairs},
                    "PRIMER_OPT_SIZE": 20, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
                    "PRIMER_OPT_TM": ${primer_opt_tm},
                    "PRIMER_MIN_TM": ${primer_min_tm},
                    "PRIMER_MAX_TM": ${primer_max_tm},
                    "PRIMER_MIN_GC": ${primer_min_gc},
                    "PRIMER_MAX_GC": ${primer_max_gc},
                    "PRIMER_PRODUCT_SIZE_RANGE": [[${product_size_min}, ${product_size_max}]],
                },
            )
        except Exception as exc:
            print(f"Warning: Primer3 failed for {rid}: {exc}")
            continue

        for i in range(result.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            primers.append(PrimerPair(
                region_id=rid, contig=contig, region_start=start, region_end=end,
                fwd_seq=result[f"PRIMER_LEFT_{i}_SEQUENCE"],
                rev_seq=result[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                fwd_tm=result[f"PRIMER_LEFT_{i}_TM"],
                rev_tm=result[f"PRIMER_RIGHT_{i}_TM"],
                fwd_gc=result[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                rev_gc=result[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                product_size=result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                penalty=result[f"PRIMER_PAIR_{i}_PENALTY"],
            ))

    print(f"Designed {len(primers)} primer pairs")

    # ── Write outputs ────────────────────────────────────────────────────────
    fields = list(PrimerPair.__dataclass_fields__.keys())
    with open(f"{prefix}_primers.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\\t")
        w.writeheader()
        for pp in primers:
            w.writerow(asdict(pp))

    with open(f"{prefix}_primers.json", "w") as fh:
        json.dump([asdict(pp) for pp in primers], fh, indent=2)

    with open(f"{prefix}_primers.fasta", "w") as fh:
        for pp in primers:
            fh.write(f">{pp.region_id}_fwd\\n{pp.fwd_seq}\\n")
            fh.write(f">{pp.region_id}_rev\\n{pp.rev_seq}\\n")

    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n')
        fh.write(f"    primer3-py: {primer3.__version__}\\n")
        import Bio
        fh.write(f"    biopython: {Bio.__version__}\\n")
    """
}
