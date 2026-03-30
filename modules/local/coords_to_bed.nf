process COORDS_TO_BED {
    tag "${meta.id}"
    label 'process_single'

    // Pure awk — any base container works
    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(coords)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix             = task.ext.prefix ?: "${meta.id}"
    def min_identity   = task.ext.min_identity   ?: 80.0
    def min_aln_length = task.ext.min_aln_length ?: 100
    """
    # The nf-core nucmer module produces show-coords pipe-delimited output:
    #
    #     [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
    #     =====================================================================================
    #      259    21548  |        1    21291  |    21290    21291  |    99.99  | REF_TAG\\tQRY_TAG
    #
    # Extract: ref_contig (S1-based coords), convert to 0-based BED,
    # filter by %identity and alignment length.

    awk -F'|' '
        # Skip header lines (before and including the === separator)
        /^=/ { data = 1; next }
        !data { next }
        NF >= 4 {
            # Field 1: "   S1     E1  "
            # Field 2: "   S2     E2  "
            # Field 3: "  LEN1   LEN2  "
            # Field 4: "  %IDY  "
            # Field 5: " REF_TAG\\tQRY_TAG"

            n1 = split(\$1, a, " +")
            n3 = split(\$3, c, " +")
            n4 = split(\$4, d, " +")

            s1 = a[2] + 0
            e1 = a[3] + 0
            len1 = c[2] + 0
            idy = d[2] + 0

            # Tags field: " REF_TAG\\tQRY_TAG" or " REF_TAG  QRY_TAG"
            n5 = split(\$5, tags, /[\\t ]+/)
            ref_tag = tags[2]

            if (idy >= ${min_identity} && len1 >= ${min_aln_length} && ref_tag != "") {
                start = s1; end = e1
                if (start > end) { tmp = start; start = end; end = tmp }
                print ref_tag "\\t" start - 1 "\\t" end
            }
        }
    ' ${coords} | sort -k1,1 -k2,2n > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -W version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
