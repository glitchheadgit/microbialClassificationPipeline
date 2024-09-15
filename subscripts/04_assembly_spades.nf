include {spades} from "./process.nf"


workflow assembly {
    take:
    index
    reads

    main:
    spades(index, reads)
}

workflow {
    c = 1
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .map { tuple(c++, it[1]) }
        .multiMap {
            index: it[0]
            reads: it[1]
        }
        .set { reads_index_ch }
    assembly( reads_index_ch.index.collect(), reads_index_ch.reads.collect(flat: false) )
}
