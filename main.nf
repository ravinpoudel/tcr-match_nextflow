#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


//include
include {

prepare_fasta;
map_TCRseq_to_TCRdb;
generate_kmers_for_TCRMatch;
run_prodigal;
generate_kmers_from_genome_protein;
cross_genomes_TCRMatch_kmers;

} from './modules.nf'

//workfolow

workflow {

    //check outdir and build it if it does not exists
    if (!params.outdir) {
        params.outdir = "${PWD}/Result"
    }
    if(!file(params.outdir).exists()){
        file(params.outdir).mkdirs()
    }
    println(sprintf("Output dir: %s", file(params.outdir)))

      //call BGC from the reads dir
	if (params.tcrbetaseq && params.genomes_dir && params.tcrmatch_database){
        
        samples = Channel.fromPath("${params.tcrbetaseq}/*.txt")
        // samples.view()
	
		//PREPARE DATA
		prepare_fasta(samples, params.outdir)
        // prepare_fasta.out.view()

        // RUN TCRMatch
        tcrmatch_database = params.tcrmatch_database
        map_TCRseq_to_TCRdb(prepare_fasta.out, tcrmatch_database,  params.outdir)
        generate_kmers_for_TCRMatch(map_TCRseq_to_TCRdb.out, "${projectDir}/scripts/", params.outdir)


        // RUN prodigal and parse into kmers from genomes
        genomes = Channel.fromPath("${params.genomes_dir}/*.f*a")
        run_prodigal(genomes,  params.outdir)
        generate_kmers_from_genome_protein(run_prodigal.out, "${projectDir}/scripts/", params.outdir)

        zz = generate_kmers_for_TCRMatch.out.combine(generate_kmers_from_genome_protein.out)
        zz.view()
    
        // Run kmer match between the results from TCRMatch and genomes
        cross_genomes_TCRMatch_kmers(zz, "${projectDir}/scripts/", params.outdir)


    }
}

