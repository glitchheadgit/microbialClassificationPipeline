include {metawrap_binrefinement}  from "./process.nf"

workflow binrefinement_metawrap {
    take:
    metabat2_bins
    concoct_bins
    maxbin2_bins

    main:
    metawrap_binrefinement( metabat2_bins, concoct_bins, maxbin2_bins )

    emit:
    refined_bins = metawrap_binrefinement.out.refined_bins
}

workflow {
    bins = params.bins.split(',')
    Channel
        .fromPath( bins.toList(), checkIfExists: true )
        .branch {
            metabat2: it ==~ /.*metabat2_bins.*/
            concoct: it ==~ /.*concoct_bins.*/
            maxbin2: it ==~ /.*maxbin2_bins.*/
        }
        .set { bins_ch }
    binrefinement_metawrap( bins_ch.metabat2.collect(), bins_ch.concoct.collect(), bins_ch.maxbin2.collect() )
}
