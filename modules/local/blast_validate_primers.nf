process BLAST_VALIDATE_PRIMERS {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/blast:2.16.0--hc155240_2'

    input:
    tuple val(meta), path(primers_fasta), path(primers_tsv)
    path blast_db_dir

    output:
    tuple val(meta), path("${prefix}_primers_validated.tsv"), emit: validated_tsv
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Find the BLAST db prefix from the .nal or .ndb file in the directory
    DB_PREFIX=\$(ls ${blast_db_dir}/*.nal 2>/dev/null | head -1 | sed 's/\\.nal\$//')
    if [ -z "\$DB_PREFIX" ]; then
        DB_PREFIX=\$(ls ${blast_db_dir}/*.ndb 2>/dev/null | head -1 | sed 's/\\.ndb\$//')
    fi

    blastn \\
        -query ${primers_fasta} \\
        -db \$DB_PREFIX \\
        -out blast_hits.tsv \\
        -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue' \\
        -evalue 10 \\
        -word_size 7 \\
        -task blastn-short \\
        -num_threads ${task.cpus} \\
        -max_target_seqs 10

    # Flag primers with >=90% identity over >=90% of primer length
    awk -F'\\t' '\$3 >= 90.0 && \$4 >= 0.9 * \$5 { print \$1 }' blast_hits.tsv | \\
        sort -u > bad_primers.txt

    # Filter the TSV
    python3 -c "
import csv, sys
bad_file = 'bad_primers.txt'
bad = set()
with open(bad_file) as f:
    content = f.read().strip()
    if content:
        bad = set(content.split('\\n'))

with open('${primers_tsv}') as fin, open('${prefix}_primers_validated.tsv', 'w', newline='') as fout:
    reader = csv.DictReader(fin, delimiter='\\t')
    writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter='\\t')
    writer.writeheader()
    kept = 0
    for row in reader:
        rid = row['region_id']
        if f'{rid}_fwd' not in bad and f'{rid}_rev' not in bad:
            writer.writerow(row)
            kept += 1
    print(f'BLAST validation: kept {kept} primer pairs', file=sys.stderr)
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -1 | sed 's/blastn: //')
    END_VERSIONS
    """
}
