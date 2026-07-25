cmake_minimum_required(VERSION 3.16)

# ── test_ngpmcmc ───────────────────────────────────────────────────────────────
add_executable(test_ngpmcmc test_ngpmcmc.c)
target_link_libraries(test_ngpmcmc PRIVATE lafp_shared)

# test_ngpmcmc requires y.txt and tobs.txt data files in the working directory;
# mark as a manual test for now (WILL_FAIL would be misleading without data).
add_test(
    NAME    test_ngpmcmc
    COMMAND test_ngpmcmc
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(test_ngpmcmc PROPERTIES
    WILL_FAIL FALSE
    LABELS    "integration"
)

# ── test_heinkel ───────────────────────────────────────────────────────────────
add_executable(test_heinkel test_heinkel.c)
target_link_libraries(test_heinkel PRIVATE lafp_shared)

add_test(
    NAME    test_heinkel
    COMMAND test_heinkel
)
set_tests_properties(test_heinkel PROPERTIES LABELS "unit")
