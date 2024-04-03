include { TRIM; QC; INDEX; HISAT; SAMTOBAM; SORT; FEATURECOUNTS; SRADOWNLOAD; OPTIMIZEDMAPPER } from './modules.nf'

params.species = ''
params.annot_format = 'gff3'
params.accession_list = ''
params.optimized = false

accession_ch = Channel
                    .fromPath(params.accession_list)
                    .splitCsv(header:false)

if (params.annot_format == 'gff3')
    annotation_ch = Channel.fromPath("$params.species*.gff3")
else
    annotation_ch = Channel.fromPath("$params.species*.gtf")


workflow  {
    
    if (params.optimized) {
        genome_ch = Channel.fromPath("$params.species*.fa")
        index_ch = INDEX(genome_ch)
        optimized_ch = OPTIMIZEDMAPPER(accession_ch, params.annot_format, annotation_ch.first(), index_ch.first())
        }
    else {
        genome_ch = Channel.fromPath("$params.species*.fa")
        index_ch = INDEX(genome_ch)
        sra_ch = SRADOWNLOAD(accession_ch)

        trimmomatic_ch = TRIM(sra_ch)
        trimmomatic_ch.view()

        fastqc_ch = QC(trimmomatic_ch)
        fastqc_ch.view()

        hisat2_ch = HISAT(trimmomatic_ch, index_ch.first())
        hisat2_ch.view()

        samtools_ch = SAMTOBAM(hisat2_ch)
        samtools_ch.view()

        sort_ch = SORT(samtools_ch)
        sort_ch.view()

        featurecounts_ch = FEATURECOUNTS(sort_ch, params.annot_format, annotation_ch.first())
        featurecounts_ch.view()
        }
}

