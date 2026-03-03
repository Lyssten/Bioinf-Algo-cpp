cmake_minimum_required(VERSION 3.10)
project(BioinfAlgo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Lecture 01: String algorithms
add_subdirectory(lectures/lecture_01/pattern_count)
add_subdirectory(lectures/lecture_01/reverse_complement)
add_subdirectory(lectures/lecture_01/frequent_words)
add_subdirectory(lectures/lecture_01/pattern_matching)
add_subdirectory(lectures/lecture_01/clump_finding)
add_subdirectory(lectures/lecture_01/minimum_skew)
add_subdirectory(lectures/lecture_01/hamming_distance)
add_subdirectory(lectures/lecture_01/frequent_words_mismatches)

# Lecture 02: Encoding and approximate matching
add_subdirectory(lectures/lecture_02/pattern_to_number)
add_subdirectory(lectures/lecture_02/number_to_pattern)
add_subdirectory(lectures/lecture_02/approximate_pattern_matching)
add_subdirectory(lectures/lecture_02/frequency_array)
add_subdirectory(lectures/lecture_02/d_neighborhood)
add_subdirectory(lectures/lecture_02/frequent_words_mismatches_complements)

# Lecture 03: Motif finding
add_subdirectory(lectures/lecture_03/motif_enumeration)
add_subdirectory(lectures/lecture_03/median_string)
add_subdirectory(lectures/lecture_03/profile_most_probable)
add_subdirectory(lectures/lecture_03/distance_pattern_strings)
add_subdirectory(lectures/lecture_03/greedy_motif_search)
add_subdirectory(lectures/lecture_03/randomized_motif_search)
add_subdirectory(lectures/lecture_03/gibbs_sampler)
