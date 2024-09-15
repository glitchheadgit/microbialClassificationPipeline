include {fastp}  from "./process.nf"
include {multiqc}  from "./process.nf"
include {fastqc as fqc_before; fastqc as fqc_after}  from "./process.nf"


workflow preprocess {
    take:
    reads_ch 
    
    main:
    fqc_before( reads_ch, 'qc_before' )
    fastp( reads_ch )
    fqc_after( fastp.out.trimmed, 'qc_after' )

    //merge reports from multiple processes
    reports = fqc_after.out.reports.concat(fqc_before.out.reports, fastp.out.reports).collect()
    multiqc( reports, 'multiqc' )

    emit:
    trimmed = fastp.out.trimmed
}


workflow {
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    preprocess( reads_ch )
}
