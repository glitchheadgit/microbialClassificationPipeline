include {bwa as map_human; bwa as map_host}  from "./process.nf"

workflow mapping {
    take:
    trimmed_reads_ch 
    ref_list

    main:
    if ( ref_list.size() == 2 ) {
        map_human( trimmed_reads_ch, ref_list[0], params.human_dir )
        map_host( map_human.out.unmapped_tuple, ref_list[1], params.host_dir )
    } else if ( ref_list.size() == 1 ) {
        map_host( trimmed_reads_ch, ref_list[0], params.host_dir )
    } else {
        println "Wrong number of references, specify it only for a host and a human(optional)."
    }
    emit:
    unmapped_tuple = map_host.out.unmapped_tuple
    
}


workflow {
    ref_list = params.ref.split(/,/)
    ref_list = ref_list.collect{ new File(it).getCanonicalPath() }
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }

    mapping( reads_ch, ref_list )
}
