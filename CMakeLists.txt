cmake_minimum_required(VERSION 3.10)
project(BioinfAlgo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

function(configure_task_output target_name relative_dir)
    set(output_dir "${PROJECT_SOURCE_DIR}/${relative_dir}")

    set_target_properties(${target_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${output_dir}"
    )

    foreach(config Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${config}" config_upper)
        set_target_properties(${target_name} PROPERTIES
            "RUNTIME_OUTPUT_DIRECTORY_${config_upper}" "${output_dir}"
        )
    endforeach()
endfunction()

function(add_bioinf_task relative_dir target_name)
    add_subdirectory(${relative_dir})
    configure_task_output(${target_name} ${relative_dir})
endfunction()

# Lecture 01: String algorithms
add_bioinf_task(lectures/lecture_01/pattern_count pattern_count)
add_bioinf_task(lectures/lecture_01/reverse_complement reverse_complement)
add_bioinf_task(lectures/lecture_01/frequent_words frequent_words)
add_bioinf_task(lectures/lecture_01/pattern_matching pattern_matching)
add_bioinf_task(lectures/lecture_01/clump_finding clump_finding)
add_bioinf_task(lectures/lecture_01/minimum_skew minimum_skew)
add_bioinf_task(lectures/lecture_01/hamming_distance hamming_distance)
add_bioinf_task(lectures/lecture_01/frequent_words_mismatches frequent_words_mismatches)

# Lecture 02: Encoding and approximate matching
add_bioinf_task(lectures/lecture_02/pattern_to_number pattern_to_number)
add_bioinf_task(lectures/lecture_02/number_to_pattern number_to_pattern)
add_bioinf_task(lectures/lecture_02/approximate_pattern_matching approximate_pattern_matching)
add_bioinf_task(lectures/lecture_02/frequency_array frequency_array)
add_bioinf_task(lectures/lecture_02/d_neighborhood d_neighborhood)
add_bioinf_task(lectures/lecture_02/frequent_words_mismatches_complements frequent_words_mismatches_complements)

# Lecture 03: Motif finding
add_bioinf_task(lectures/lecture_03/motif_enumeration motif_enumeration)
add_bioinf_task(lectures/lecture_03/median_string median_string)
add_bioinf_task(lectures/lecture_03/profile_most_probable profile_most_probable)
add_bioinf_task(lectures/lecture_03/distance_pattern_strings distance_pattern_strings)
add_bioinf_task(lectures/lecture_03/greedy_motif_search greedy_motif_search)
add_bioinf_task(lectures/lecture_03/randomized_motif_search randomized_motif_search)
add_bioinf_task(lectures/lecture_03/gibbs_sampler gibbs_sampler)
