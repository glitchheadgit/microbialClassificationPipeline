include {preprocess} from "./subscripts/01_preprocess.nf"
include {mapping} from "./subscripts/02_mapping.nf"
include {classify_reads} from "./subscripts/03_classify_reads.nf"
include {assembly_megahit} from "./subscripts/04_assembly.nf"
include {binning_metawrap} from "./subscripts/05_binning.nf"
include {binrefinement_metawrap} from "./subscripts/06_binrefinement.nf"
include {gtdbtk_classify} from "./subscripts/07_classify_bins.nf"
include {gunzip} from "./subscripts/process.nf"


workflow {
    //preprocessing
    Channel
        .fromFilePairs( params.reads, checkIfExists: true )
        .set { reads_ch }
    preprocess( reads_ch )

    //mapping
    ref_list = params.ref.split( /,/ )
    ref_list = ref_list.collect{ new File(it).getAbsolutePath() }
    mapping( preprocess.out.trimmed, ref_list )

    //reads classification
    classify_reads( mapping.out.unmapped_tuple )

    //megahit assembly
    mapping.out.unmapped_tuple
    .multiMap {
        read1: it[1][0]
        read2: it[1][1]
    }
    .set { read1_read2_ch }
    assembly_megahit( read1_read2_ch.read1.collect(), read1_read2_ch.read2.collect() )

    //metawrap binning
    gunzip(  read1_read2_ch.read1.collect(), read1_read2_ch.read2.collect()  )
    binning_metawrap( assembly_megahit.out.assembly, gunzip.out.reads1, gunzip.out.reads2 )
    //metawrap binrefinement
    binning_metawrap.out.metabat2_bins.ifEmpty( 'null1' )
    binning_metawrap.out.concoct_bins.ifEmpty( 'null2' )
    binning_metawrap.out.maxbin2_bins.ifEmpty( 'null3' )
    binrefinement_metawrap( binning_metawrap.out.metabat2_bins.collect(), binning_metawrap.out.concoct_bins.collect(), binning_metawrap.out.maxbin2_bins.collect() )
    gtdbtk_classify( binrefinement_metawrap.out.refined_bins )
}
