
process TRIM{
    tag 'Trimmomatic'
    input:
        tuple val(sample_id), val(layout), path(reads)
    
    output:
        tuple val(sample_id), val(layout), path('*trim.fastq')
    
    script:
        if (layout.toUpperCase() == "SINGLE")
            """
            trimmomatic SE -phred33 $reads ${sample_id}_trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

            """
        else 
            """
            trimmomatic PE -phred33 $reads ${sample_id}_1_trim.fastq ${sample_id}_1_trimU.fastq ${sample_id}_2_trim.fastq ${sample_id}_2_trimU.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

            """
}


process QC{
    tag 'FastQC'
    publishDir "Results/${params.species}/fastqc-reports", mode: 'copy'

    input:
        tuple val(sample_id), val(layout), path(trimmed_reads)
    
    output:
        path '*fastqc.zip'
    
    script:
        """
        fastqc $trimmed_reads
        """
}


process INDEX{
    tag 'INDEX'
    /*publishDir "Results/"+"$params.species", mode: 'copy'*/

    input:
        path genome
    
    output:
        path '*.ht2'
    
    script:
        """
        hisat2-build -p 2 -f $genome $params.species
        """
}


process HISAT{
    tag 'HISAT2'

    input:
        tuple val(sample_id), val(layout), path(reads)
        path index
    
    output:
        tuple val(sample_id), path('*.sam')

    script:
        if (layout.toUpperCase() == 'SINGLE')
            """
            hisat2 -p 5 --phred33 -x $params.species -U $reads -S ${sample_id}.sam
            """
        else
            """
            hisat2 -p 5 --phred33 -x $params.species -1 ${reads[0]} -2 ${reads[0]} -S ${sample_id}.sam
            """
}

process SAMTOBAM{
    tag 'Samtools'

    input:
        tuple val(sample_id), path(sam_file)
    
    output:
        tuple val(sample_id), path('*.unsorted.bam')
    
    script:
        """
        samtools view -b $sam_file > ${sam_file}.unsorted.bam
        """
}

process SORT{
    tag 'Sort'

    input:
        tuple val(sample_id), path(unsorted_bam)
    
    output:
        tuple val(sample_id), path('*.sorted.bam')
    
    script:
        """
        samtools sort $unsorted_bam > ${unsorted_bam}.sorted.bam
        """
}

process FEATURECOUNTS{
    tag 'FeatureCounts'
    publishDir "Results/${params.species}/read-counts", mode: 'copy'

    input:
        tuple val(sample_id), path(sorted_bam)
        val  annot_format
        path annotation
    
    output:
        path '*_counts.txt'
    
    script:
        if (annot_format == 'gff3')
            """
            featureCounts -t mRNA -g Parent -a $annotation -o ${sample_id}_counts.txt $sorted_bam
            """
        else
            """
            featureCounts -t transcript -g gene_id -a $annotation -o ${sample_id}_counts.txt $sorted_bam
            """
}

process SRADOWNLOAD{
    tag 'NCBI SRA DOWNLOAD'

    input:
        tuple val(sample_id), val(layout)
    
    output:
        tuple val(sample_id), val(layout), path('*.fastq')
    
    script:
        """
        fasterq-dump $sample_id -e 5 -p
        """
}

process OPTIMIZEDMAPPER {
    scratch 'scratch_dir'
    maxForks 3
    publishDir "Results/${params.species}/read-counts", mode: 'move', pattern: '*counts.txt'
    publishDir "Results/${params.species}/fastqc-reports", mode: 'move', pattern: '*fastqc.zip'


    input:
        tuple val(sample_id), val(layout)
        val  annot_format
        path annotation
        path hisat2_index
    
    output:
        path '*fastqc.zip'
        path '*_counts.txt'
    
    script:
        if (layout.toUpperCase() == 'SINGLE') {
            """
            fasterq-dump $sample_id -e 5 -p
            trimmomatic SE -phred33 ${sample_id}.fastq ${sample_id}_trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
            fastqc ${sample_id}_trim.fastq
            hisat2 -p 5 --phred33 -x $params.species -U ${sample_id}_trim.fastq -S ${sample_id}.sam
            samtools view -b ${sample_id}.sam > ${sample_id}.unsorted.bam
            samtools sort ${sample_id}.unsorted.bam > ${sample_id}.bam
            featureCounts -t mRNA -g Parent -a $annotation -o ${sample_id}_counts.txt ${sample_id}.bam
            """
            }
        else {
           """
            fasterq-dump $sample_id -e 5 -p
            trimmomatic PE -phred33 ${sample_id}_1.fastq ${sample_id}_2.fastq ${sample_id}_1_trim.fastq ${sample_id}_1_trimU.fastq ${sample_id}_2_trim.fastq ${sample_id}_2_trimU.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
            fastqc ${sample_id}_1_trim.fastq ${sample_id}_2_trim.fastq
            hisat2 -p 5 --phred33 -x $params.species -1 ${sample_id}_1_trim.fastq -2 ${sample_id}_2_trim.fastq -S ${sample_id}.sam
            samtools view -b ${sample_id}.sam > ${sample_id}.unsorted.bam
            samtools sort ${sample_id}.unsorted.bam > ${sample_id}.bam
            featureCounts -t mRNA -g Parent -a $annotation -o ${sample_id}_counts.txt ${sample_id}.bam
           """
            }
}
    
