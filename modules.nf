#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process prepare_fasta {

    tag "$prefix"
    errorStrategy 'retry'
    maxRetries 5

    input:
        file tcrseq_file
        path workdir
    output:
        tuple val(prefix), path("Samples/${prefix}/prepare_tcrseq")

    storeDir "${workdir}"

    script:
    prefix=tcrseq_file.name.split('.txt')[0]
    """
    mkdir -p Samples/${prefix}/prepare_tcrseq
    cp $tcrseq_file Samples/${prefix}/prepare_tcrseq/${prefix}.txt
    """

}


process map_TCRseq_to_TCRdb{
    tag "$prefix"
    errorStrategy 'retry'
    maxRetries 5

    input:
        tuple val(prefix), path(tcrseq_path)
        path iedb_data
        path workdir

    output:
        tuple  val(prefix), path("Samples/${prefix}/tcrmatch")
    
    storeDir "${workdir}"
    
    script:    
    """
        mkdir -p Samples/${prefix}/tcrmatch/tmp
        cp $tcrseq_path/*.txt  Samples/${prefix}/tcrmatch/tmp/tcrseq.txt
        /app/tcrmatch -i  Samples/${prefix}/tcrmatch/tmp/tcrseq.txt -t 32 -d $iedb_data -a  > Samples/${prefix}/tcrmatch/TCRMatch.txt

    """
}



process generate_kmers_for_TCRMatch{
    tag "${prefix}"
    errorStrategy 'retry'
    maxRetries 5

    input:
        tuple val(prefix), path(tcrmatchpath)
        path bindir
        path workdir

    output:
        tuple  val(prefix), path("Samples/${prefix}/tcrmatch_to_kmers")
    
    storeDir "${workdir}"
    
    script:
   
    """
        mkdir -p Samples/${prefix}/tcrmatch_to_kmers
        python3 $bindir/process_tcrmatch_kmers.py $tcrmatchpath/TCRMatch.txt Samples/${prefix}/tcrmatch_to_kmers/tcrmatch_8_11_mers.txt

    """
}



process run_prodigal{
    tag "${sample}"
    errorStrategy 'retry'
    maxRetries 5

    input:
        path input_genome
        path workdir

    output:
        path "Genomes/genomes_to_faa/${input_genome.baseName}.faa"
    
    storeDir "${workdir}"
    
    script:
    sample=input_genome.name.split('.fa')[0]
    """
        echo $sample
        mkdir -p Genomes/genomes_to_faa
        pyrodigal -i $input_genome -o Genomes/genomes_to_faa/${input_genome.baseName}.genes -a Genomes/genomes_to_faa/${input_genome.baseName}.faa

    """
}

process generate_kmers_from_genome_protein{
    errorStrategy 'retry'
    maxRetries 5

    input:
        file genomes_faa
        path bindir
        path workdir

    output:
        path("Genomes/genomes_to_kmers/${genomes_faa.baseName}_genomes_8_11_mers.txt")
    
    storeDir "${workdir}"
    
    script:
   
    """
        mkdir -p Genomes/genomes_to_kmers
        python3 $bindir/process_genomes.py $genomes_faa Genomes/genomes_to_kmers/${genomes_faa.baseName}_genomes_8_11_mers.txt

    """
}

    
process cross_genomes_TCRMatch_kmers{
    tag "${prefix}"
    errorStrategy 'retry'
    maxRetries 5

    input:
        tuple val(prefix), path(tcrmatch_kmers), file(genome_kmers)
        path bindir
        path workdir

    output:
        file("Genomes/cross_TCRMatch_Genomes_Epitopes/${prefix}_${genome_kmers.baseName}_tcrmatch_genomes_cross_match.txt")
        file("Genomes/cross_TCRMatch_Genomes_Epitopes/${prefix}_${genome_kmers.baseName}_tcrmatch_genomes_cross_match_summary.txt")
        file("Genomes/cross_TCRMatch_Genomes_Epitopes/${prefix}_${genome_kmers.baseName}_tcrmatch_genomes_cross_match_only_crags.txt")
        

    storeDir "${workdir}"
    
    script:
   
    """
        mkdir -p Genomes/cross_TCRMatch_Genomes_Epitopes
        Rscript $bindir/kmer_match_detailed_cl_many_kmers.R \
        $tcrmatch_kmers/tcrmatch_8_11_mers.txt \
        $genome_kmers \
        Genomes/cross_TCRMatch_Genomes_Epitopes/${prefix}_${genome_kmers.baseName}_tcrmatch_genomes_cross_match.txt \
        Genomes/cross_TCRMatch_Genomes_Epitopes/${prefix}_${genome_kmers.baseName}_tcrmatch_genomes_cross_match_summary.txt \
        2 \
        Genomes/cross_TCRMatch_Genomes_Epitopes/${prefix}_${genome_kmers.baseName}_tcrmatch_genomes_cross_match_only_crags.txt 

    """
}


