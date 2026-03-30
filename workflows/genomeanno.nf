/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genomeanno_pipeline'
include { GTDBTK_CLASSIFYWF } from '../modules/nf-core/gtdbtk/classifywf/main'
include { ABRICATE_RUN } from '../modules/nf-core/abricate/run/main'
include { ABRICATE_SUMMARY } from '../modules/nf-core/abricate/summary/main'
include { CHECKM2_PREDICT } from '../modules/nf-core/checkm2/predict/main'
include { NUCMER                 } from '../modules/nf-core/nucmer/main'
include { EXTRACT_CLOSE_REFS     } from '../modules/local/extract_close_refs'
include { COORDS_TO_BED          } from '../modules/local/coords_to_bed'
include { DESIGN_STRAIN_PRIMERS  } from '../modules/local/design_strain_primers'
include { BLAST_VALIDATE_PRIMERS } from '../modules/local/blast_validate_primers'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEANNO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Abricate
    //
    if (!params.arg_skip_abricate) {
    def abricate_dbs = params.arg_abricate_db_ids instanceof List
        ? params.arg_abricate_db_ids
        : [params.arg_abricate_db_ids]

    // Create a channel of [db_name] values
    ch_abricate_dbs = Channel.fromList(abricate_dbs)

    // Combine every sample with every DB → [meta + db_id, assembly, db]
    ch_abricate_input = ch_samplesheet
        .combine(ch_abricate_dbs)
        .map { meta, assembly, db ->
            def new_meta = meta + [db: db]
            [ new_meta, assembly ]
        }

    // Use the custom databasedir if provided
    def abricate_db_dir = params.arg_abricate_db ? file(params.arg_abricate_db) : []

    ABRICATE_RUN(
        ch_abricate_input,
        abricate_db_dir
    )
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions.first())

    // Group by DB for per-DB summary
    ABRICATE_RUN.out.report
        .map { meta, report -> [ meta.db, meta, report ] }
        .groupTuple(by: 0)
        .map { db, metas, reports ->
            [ [ id: "summary_${db}", db: db ], reports ]
        }
        | ABRICATE_SUMMARY

    ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions.first())
    }
    //
    // MODULE: GTDB-Tk (Taxonomic Classification)
    //
    if (params.gtdb_db) {
        // Create a value channel for the database so it can be reused for every sample
        ch_gtdb_db = channel.fromPath(params.gtdb_db).map{ db -> [ [id:'gtdb'], db ] }.first()

        GTDBTK_CLASSIFYWF (
            ch_samplesheet,
            ch_gtdb_db,
            false
        )       
        ch_versions      = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(GTDBTK_CLASSIFYWF.out.summary.collect{v -> v[1]})
    }

    //
    // STRAIN-SPECIFIC PRIMER DESIGN
    // Requires: GTDB-Tk with ANI screen enabled (default)
    //
    if (params.gtdb_db && !params.skip_primer_design && !params.gtdbtk_skip_aniscreen) {
        
        if (params.primer_gtdbtk_outdir) {
            ch_gtdbtk_for_primers = ch_samplesheet
                .map { meta, _assembly ->
                    def sample_dir = file("${params.primer_gtdbtk_outdir}/${meta.id}")
                    def gtdbtk_dir = sample_dir.exists() ? sample_dir : file(params.primer_gtdbtk_outdir)
                    [ meta, gtdbtk_dir ]
                }
        } else if (!params.gtdbtk_skip_aniscreen) {
            ch_gtdbtk_for_primers = GTDBTK_CLASSIFYWF.out.gtdb_outdir
        } else {
            ch_gtdbtk_for_primers = Channel.empty()
        }
        //
        // Step 1: Extract closest reference genomes from GTDB-Tk ANI output
        //
        EXTRACT_CLOSE_REFS(
            GTDBTK_CLASSIFYWF.out.gtdb_outdir,       // [meta, gtdbtk_output_dir/]
            file(params.gtdb_db, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(EXTRACT_CLOSE_REFS.out.versions.first())

        //
        // Step 2: Align each reference genome against the target assembly (nucmer)
        //
        ch_nucmer_input = ch_samplesheet
            .join(EXTRACT_CLOSE_REFS.out.refs)        // [meta, assembly, [ref1.fna, ref2.fna, ...]]
            .flatMap { meta, assembly, refs ->
                def asm = assembly instanceof List ? assembly[0] : assembly
                refs.collect { ref ->
                    def nucmer_meta = meta + [ref_name: ref.baseName]
                    [ nucmer_meta, asm, ref ]
                }
            }

        NUCMER(ch_nucmer_input)
        ch_versions = ch_versions.mix(NUCMER.out.versions.first())

        //
        // Step 3: Convert coords to BED, then collect all BEDs per sample
        //
        COORDS_TO_BED(NUCMER.out.coords)
        ch_versions = ch_versions.mix(COORDS_TO_BED.out.versions.first())

        // Group BED files back by original sample ID (strip ref_name from meta)
        ch_beds_per_sample = COORDS_TO_BED.out.bed
            .map { meta, bed ->
                def sample_meta = meta.subMap('id') + meta.findAll { k, _v -> k != 'ref_name' }
                // Use just the sample id for grouping
                [ meta.id, meta, bed ]
            }
            .groupTuple(by: 0)
            .map { sample_id, metas, beds ->
                // Recover the original meta (without ref_name)
                def original_meta = metas[0].findAll { k, _v -> k != 'ref_name' }
                [ original_meta, beds ]
            }

        //
        // Step 4: Find unique regions and design primers
        //
        ch_primer_input = ch_samplesheet
            .map { meta, assembly ->
                def asm = assembly instanceof List ? assembly[0] : assembly
                [ meta.id, meta, asm ]
            }
            .join(
                ch_beds_per_sample.map { meta, beds -> [ meta.id, beds ] },
                by: 0
            )
            .map { sample_id, meta, assembly, beds ->
                [ meta, assembly, beds ]
            }

        DESIGN_STRAIN_PRIMERS(ch_primer_input)
        ch_versions = ch_versions.mix(DESIGN_STRAIN_PRIMERS.out.versions.first())

        //
        // Step 5 (optional): BLAST validation
        //
        if (params.primer_blast_db) {
            ch_blast_input = DESIGN_STRAIN_PRIMERS.out.primers_fasta
                .join(DESIGN_STRAIN_PRIMERS.out.primers_tsv)
                .map { meta, fasta, tsv -> [ meta, fasta, tsv ] }

            BLAST_VALIDATE_PRIMERS(
                ch_blast_input,
                file(params.primer_blast_db, checkIfExists: true)
            )
            ch_versions = ch_versions.mix(BLAST_VALIDATE_PRIMERS.out.versions.first())
        }
    }


    if (params.checkm2_db) {
        ch_checkm2_db = [[:], file(params.checkm2_db, checkIfExists: true)]

        CHECKM2_PREDICT(ch_samplesheet, ch_checkm2_db)
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(
            CHECKM2_PREDICT.out.checkm2_tsv.map { _meta, summary -> summary }.flatten()
        )
    }
    
    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'genomeanno_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
