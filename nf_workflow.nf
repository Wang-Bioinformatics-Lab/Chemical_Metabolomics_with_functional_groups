#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.inputlibraries = "data/libraries"
params.inputspectra = "data/spectra"

// Parameters
params.searchtool = "gnps" // blink, gnps

params.topk = 500

params.fragment_tolerance = 0.5
params.pm_tolerance = 2.0

params.library_min_cosine = 0.7
params.library_min_matched_peaks = 6

params.merge_batch_size = 1000 //Not a UI parameter

// Filtering structures
params.filtertostructures = "1" // 1 means we filter to only hits with structures

//TODO: Implement This
params.filter_precursor = 1
params.filter_window = 1

//TODO: Implement This
params.analog_search = "1"
params.analog_max_shift = 1999

// Blink Parameters
params.blink_ionization = "positive"
params.blink_minpredict = 0.01

include {searchDataGNPS} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {searchDataBlink} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {formatBlinkResults} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {chunkResults} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {mergeResults} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {librarygetGNPSAnnotations} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {filtertop1Annotations} from './bin/LibrarySearch_Workflow/nf_workflow.nf'
include {summaryLibrary} from './bin/LibrarySearch_Workflow/nf_workflow.nf'

params.with_library_search = "0"
params.and_order_strategy = "0"
params.remove_knowns = "1"
params.library_search_result = "data/merged_results_with_gnps.tsv"
params.filter_threshold = 0.5
params.just_show_results = "0"
params.just_results_file = ""

TOOL_FOLDER = "$baseDir/bin"

process Extract_FGs {
    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library 
    path concatfgs

    output:
    file 'temp/*.csv'

    """
    mkdir temp
    python $TOOL_FOLDER/extract_fgs.py $library temp
    """
}

process Merge_FGs {
    publishDir "./nf_output", mode: 'copy' 

    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path concatfgs 

    output:
    file 'fgs.csv'
    file 'images/*'

    """
    mkdir images
    python $TOOL_FOLDER/merge_fgs.py $concatfgs fgs.csv --images_dir images
    """
}

process Extract_Knowns {
    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library
    file passedfgs

    output:
    file 'temp/*.csv'

    """
    mkdir temp
    python $TOOL_FOLDER/extract_knowns.py $library $passedfgs temp
    """
}

process Extract_Reactants {
    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library
    file passedfgs

    output:
    file 'temp/*.csv'

    """
    mkdir temp
    python $TOOL_FOLDER/extract_reactants.py $library $passedfgs temp
    """
}

process Remove_Knowns {
    publishDir "./nf_output", mode: 'copy' 

    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library_search_result
    file knowns

    output:
    file 'merged_results_with_gnps_removed_knowns.tsv'

    """
    python $TOOL_FOLDER/remove_knowns.py $library_search_result $knowns merged_results_with_gnps_removed_knowns.tsv --topk 30
    """
}

process Calculate_Adjusted_Ranking_and {
    publishDir "./nf_output", mode: 'copy' 

    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library_search_result
    file reactants
    file fgs

    output:
    file 'merged_results_with_gnps_adjusted_ranking.tsv'
    file 'merged_results_with_gnps_adjusted_ranking_reactants_only.tsv'

    """
    python $TOOL_FOLDER/adjusted_ranking.py $library_search_result $reactants $fgs merged_results_with_gnps_adjusted_ranking.tsv \
     --reactants_only_path merged_results_with_gnps_adjusted_ranking_reactants_only.tsv  \
        --and_order_strategy
    """
}

process Calculate_Adjusted_Ranking_or {
    publishDir "./nf_output", mode: 'copy' 

    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library_search_result
    file reactants
    file fgs

    output:
    file 'merged_results_with_gnps_adjusted_ranking.tsv'
    file 'merged_results_with_gnps_adjusted_ranking_reactants_only.tsv'

    """
    python $TOOL_FOLDER/adjusted_ranking.py $library_search_result $reactants $fgs merged_results_with_gnps_adjusted_ranking.tsv \
     --reactants_only_path merged_results_with_gnps_adjusted_ranking_reactants_only.tsv
    """
}

process Calculate_Analysis_results {
    publishDir "./nf_output", mode: 'copy' 

    conda "$TOOL_FOLDER/conda_env.yml"
    // conda "/home/user/mambaforge/envs/test_gwen"

    input:
    path library_search_result
    file knowns

    output:
    file 'analysis_results/enriched_search_results.csv'
    file 'analysis_results/analysis_results.csv'
    file 'analysis_results/improved_results.csv'
    file 'analysis_results/improved_results_best.csv'
    file 'images/*'

    // if params.and_order_strategy == "1": then pass the argument --and_order_strategy
    """
    mkdir analysis_results
    mkdir images
    python $TOOL_FOLDER/get_analysis_result.py $library_search_result $knowns analysis_results images --filter_threshold $params.filter_threshold
    """
}

process Show_Results{
    publishDir "./", mode: 'copy'

    input:
    file input 

    output:
    file 'nf_output/*'

    """
    tar -xvf $input
    """
}

workflow Adjust_and_Analysis{
    take:
    annotation_results_ch_value_channel

    main:
    spectra_2 = Channel.fromPath(params.inputspectra + "/**", relative: false)
    extracted_fgs = Extract_FGs(spectra_2, annotation_results_ch_value_channel).collectFile(name: './nf_output/rep_fgs', newLine: false, keepHeader: true, skip: 1)
    (merged_fgs, fg_images) = Merge_FGs(extracted_fgs)
    merged_fgs = merged_fgs.collect()

    extracted_knowns = Extract_Knowns(spectra_2, merged_fgs).collectFile(name: './nf_output/knowns.csv', newLine: false, keepHeader: true, skip: 1).collect()

    extracted_reactants = Extract_Reactants(spectra_2, merged_fgs).collectFile(name: './nf_output/reactants.csv', newLine: false, keepHeader: true, skip: 1).collect()


    if (params.remove_knowns == "1"){
        removed_knowns = Remove_Knowns(annotation_results_ch_value_channel, extracted_knowns).collect()
    }
    else{
        removed_knowns = annotation_results_ch_value_channel
    }

    // Calculate_Adjusted_Ranking
    if (params.and_order_strategy == "1"){
        (adjusted_ranking, reactants_only) = Calculate_Adjusted_Ranking_and(removed_knowns, extracted_reactants, merged_fgs)
    }
    else{
        (adjusted_ranking, reactants_only) = Calculate_Adjusted_Ranking_or(removed_knowns, extracted_reactants, merged_fgs)
    }

    // Calculate_Analysis_results
    (enriched_search_results, analysis_results, improved_results, improved_results_best, images) = Calculate_Analysis_results(adjusted_ranking, extracted_knowns)
}

workflow Library_Search{
    main:
    libraries_ch = Channel.fromPath(params.inputlibraries + "/*.mgf" )
    spectra = Channel.fromPath(params.inputspectra + "/**", relative: true)

    // Lets create a summary for the library files
    library_summary_ch = summaryLibrary(libraries_ch)

    // Merging all these tsv files from library_summary_ch within nextflow
    library_summary_merged_ch = library_summary_ch.collectFile(name: "./nf_output/library_summary.tsv", keepHeader: true)
    
    if(params.searchtool == "gnps"){
        // Perform cartesian product producing all combinations of library, spectra
        inputs = libraries_ch.combine(spectra)

        // For each path, add the path as a string for file naming. Result is [library_file, spectrum_file, spectrum_path_as_str]
        // Must add the prepend manually since relative does not include the glob.
        inputs = inputs.map { it -> [it[0], file(params.inputspectra + '/' + it[1]), it[1].toString().replaceAll("/","_"), it[1]] }

        (search_results) = searchDataGNPS(inputs)

        chunked_results = chunkResults(search_results.buffer(size: params.merge_batch_size, remainder: true))
    
        // Collect all the batched results and merge them at the end
        merged_results = mergeResults(chunked_results.collect())
    }
    else if (params.searchtool == "blink"){
        // Must add the prepend manually since relative does not inlcude the glob.
        spectra = spectra.map { it -> file(params.inputspectra + '/' + it) }
        search_results = searchDataBlink(libraries_ch, spectra)

        formatted_results = formatBlinkResults(search_results)

        merged_results = mergeResults(formatted_results.collect())
    }
    annotation_results_ch = librarygetGNPSAnnotations(merged_results, library_summary_merged_ch)
    annotation_results_ch_value_channel = annotation_results_ch.collect()

    emit:
    annotation_results_ch_value_channel
}

// TODO: This main will define the workflow that can then be imported, not really used right now, but can be
workflow {
    if (params.just_show_results == "1" || params.just_show_results == 1){
        just_results_file = Channel.fromPath(params.just_results_file).collect()
        just_results_file.view()
        just_results_file.subscribe { println it }
        Show_Results(just_results_file)
    }
    else{
        if (params.with_library_search == "1" || params.with_library_search == 1){
            annotation_results_ch_value_channel = Library_Search()
        }
        else{
            annotation_results_ch_value_channel = Channel.fromPath(params.library_search_result, type: 'file').collectFile(name: "./nf_output/merged_results_with_gnps.tsv")
        }

        Adjust_and_Analysis(annotation_results_ch_value_channel)
    }
}