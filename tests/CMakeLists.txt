cmake_minimum_required(VERSION 3.16)

# ── test_ngpmcmc (Integration Test) ──────────────────────────────────────────
add_executable(test_ngpmcmc test_ngpmcmc.c)
target_link_libraries(test_ngpmcmc PRIVATE lafp_shared)

add_test(
    NAME    test_ngpmcmc
    COMMAND test_ngpmcmc
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(test_ngpmcmc PROPERTIES
    WILL_FAIL FALSE
    LABELS    "integration"
)

# ── test_heinkel (Unit Test) ──────────────────────────────────────────────────
add_executable(test_heinkel test_heinkel.c)
target_link_libraries(test_heinkel PRIVATE lafp_shared)

add_test(
    NAME    test_heinkel
    COMMAND test_heinkel
)
set_tests_properties(test_heinkel PROPERTIES LABELS "unit")

# ── test_kalman (Unit Test) ───────────────────────────────────────────────────
add_executable(test_kalman test_kalman.c)
target_link_libraries(test_kalman PRIVATE lafp_shared)

add_test(
    NAME    test_kalman
    COMMAND test_kalman
)
set_tests_properties(test_kalman PROPERTIES LABELS "unit")

# ── test_matrix_utils (Unit Test) ─────────────────────────────────────────────
add_executable(test_matrix_utils test_matrix_utils.c)
target_link_libraries(test_matrix_utils PRIVATE lafp_shared)

add_test(
    NAME    test_matrix_utils
    COMMAND test_matrix_utils
)
set_tests_properties(test_matrix_utils PROPERTIES LABELS "unit")

# ── macOS ad-hoc code signing for CTest targets ──────────────────────────────
if(APPLE)
    add_custom_command(TARGET test_ngpmcmc POST_BUILD
        COMMAND codesign --force --sign - $<TARGET_FILE:test_ngpmcmc>
    )
    add_custom_command(TARGET test_heinkel POST_BUILD
        COMMAND codesign --force --sign - $<TARGET_FILE:test_heinkel>
    )
    add_custom_command(TARGET test_kalman POST_BUILD
        COMMAND codesign --force --sign - $<TARGET_FILE:test_kalman>
    )
    add_custom_command(TARGET test_matrix_utils POST_BUILD
        COMMAND codesign --force --sign - $<TARGET_FILE:test_matrix_utils>
    )
endif()
