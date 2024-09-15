include {metawrapbin}  from "./process.nf"


workflow binning_metawrap {
    take:
    assembly
    reads1
    reads2

    main:
    metawrapbin( assembly, reads1, reads2 )

    emit:
    metabat2_bins = metawrapbin.out.metabat2_bins
    concoct_bins = metawrapbin.out.concoct_bins
    maxbin2_bins = metawrapbin.out.maxbin2_bins
}

workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .multiMap {
            read1: it[1][0]
            read2: it[1][1]
        }
        .set { reads_ch }
    binning_metawrap( params.assembly, reads_ch.read1.collect(), reads_ch.read2.collect() )
}
