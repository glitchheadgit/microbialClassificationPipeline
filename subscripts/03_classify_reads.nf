include {kraken2}  from "./process.nf"
include {bracken}  from "./process.nf"
include {merge_kr2_br}  from "./process.nf"
include {metaphlan}  from "./process.nf"


workflow classify_reads {
    take:
    reads

    main:
    if (params.kraken2) {
        kraken2( reads )
        bracken( kraken2.out.stats )
        merge_kr2_br( kraken2.out.stats, bracken.out.stats )
    } 
    if (params.metaphlan) {
        metaphlan( reads )
    }
}

workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    classify_reads( reads_ch )
}
