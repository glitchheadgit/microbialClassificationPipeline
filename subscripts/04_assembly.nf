include {megahit}  from "./process.nf"


workflow assembly_megahit {
    take:
    reads1
    reads2

    main:
    megahit( reads1, reads2 )
    
    emit:
    assembly = megahit.out.assembly
}

workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .multiMap {
            read1: it[1][0]
            read2: it[1][1]
        }
        .set { reads_ch }
    assembly_megahit( reads_ch.read1.collect(), reads_ch.read2.collect() )
}
