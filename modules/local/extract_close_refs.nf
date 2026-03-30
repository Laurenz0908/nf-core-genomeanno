process EXTRACT_CLOSE_REFS {
    tag "${meta.id}"
    label 'process_single'

    container 'community.wave.seqera.io/library/pip_biopython_primer3-py:e339c879b0de559a'

    input:
    tuple val(meta), path(gtdbtk_outdir)
    path gtdb_db

    output:
    tuple val(meta), path("refs/*.fna"), emit: refs
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def max_refs  = task.ext.max_refs ?: 20
    def min_ani   = task.ext.min_ani ?: 80.0
    """
    #!/usr/bin/env python3
    import csv
    import gzip
    import shutil
    import sys
    from pathlib import Path

    gtdbtk_dir = Path("${gtdbtk_outdir}")
    gtdb_db    = Path("${gtdb_db}")
    max_refs   = ${max_refs}
    min_ani    = ${min_ani}

    outdir = Path("refs")
    outdir.mkdir(exist_ok=True)

    # Find the ani_summary.tsv — it's in classify/ani_screen/
    ani_files = list(gtdbtk_dir.glob("classify/ani_screen/*.ani_summary.tsv"))
    if not ani_files:
        # Fallback: try the regular summary if ani_screen wasn't run
        ani_files = list(gtdbtk_dir.glob("classify/*.summary.tsv"))
    if not ani_files:
        sys.exit("ERROR: No GTDB-Tk summary/ani_summary file found")

    ani_file = ani_files[0]
    print(f"Parsing: {ani_file}")

    # Find genome_paths.tsv in the GTDB database (inside the skani dir)
    genome_paths_file = None
    for candidate in [
        gtdb_db / "skani" / "genome_paths.tsv",
        # Also check if gtdb_db points directly to the release directory
        *list(gtdb_db.glob("**/skani/genome_paths.tsv")),
    ]:
        if candidate.exists():
            genome_paths_file = candidate
            break

    if not genome_paths_file:
        sys.exit("ERROR: genome_paths.tsv not found in GTDB database")

    # Build lookup: genome_id -> relative path
    genome_paths = {}
    with open(genome_paths_file) as fh:
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) >= 2:
                fname = parts[0]                          # e.g. GCA_000008085.1_genomic.fna.gz
                reldir = parts[1]                         # e.g. database/GCA/000/008/085/
                genome_id = fname.replace("_genomic.fna.gz", "")
                genome_paths[genome_id] = genome_paths_file.parent / reldir / fname

    # Parse ani_summary.tsv
    ref_genomes = []

    with open(ani_file) as fh:
        reader = csv.DictReader(fh, delimiter="\\t")
        for row in reader:
            # Primary reference
            ref_id = row.get("reference_genome", "")
            if ref_id:
                ref_genomes.append((ref_id, 100.0))  # primary match

            # Other related references from the last column
            other_col = row.get(
                "other_related_references(genome_id,species_name,radius,ANI,AF)", ""
            )
            if other_col and other_col.strip():
                for entry in other_col.split("; "):
                    parts = [p.strip() for p in entry.split(",")]
                    if len(parts) >= 5:
                        gid = parts[0].strip()
                        try:
                            ani = float(parts[3].strip())
                        except ValueError:
                            continue
                        if ani >= min_ani:
                            ref_genomes.append((gid, ani))

    # Sort by ANI descending, take top N
    ref_genomes.sort(key=lambda x: x[1], reverse=True)
    ref_genomes = ref_genomes[:max_refs]

    print(f"Selected {len(ref_genomes)} reference genomes (ANI >= {min_ani})")

    # Copy and decompress reference genomes
    copied = 0
    for gid, ani in ref_genomes:
        fasta_path = genome_paths.get(gid)
        if fasta_path is None:
            # Try with RS_/GB_ prefix stripped or added
            for alt in [f"GCF{gid[3:]}", f"GCA{gid[3:]}", gid]:
                if alt in genome_paths:
                    fasta_path = genome_paths[alt]
                    break
        if fasta_path is None:
            print(f"WARNING: genome {gid} not found in genome_paths.tsv, skipping")
            continue

        if not fasta_path.exists():
            print(f"WARNING: file {fasta_path} does not exist, skipping")
            continue

        out_fna = outdir / f"{gid}_genomic.fna"
        if str(fasta_path).endswith(".gz"):
            with gzip.open(fasta_path, "rb") as fin, open(out_fna, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            shutil.copy2(fasta_path, out_fna)
        copied += 1
        print(f"  {gid} (ANI={ani:.2f})")

    if copied == 0:
        sys.exit("ERROR: No reference genomes could be extracted")

    print(f"Extracted {copied} reference genomes to {outdir}/")

    # versions
    with open("versions.yml", "w") as fh:
        fh.write('"EXTRACT_CLOSE_REFS":\\n')
        fh.write("    python: \$(python3 --version | sed 's/Python //')\\n")
    """
}
