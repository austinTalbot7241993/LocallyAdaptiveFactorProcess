cmake_minimum_required(VERSION 3.16)

# ── run_ngpmcmc ────────────────────────────────────────────────────────────────
# Usage: run_ngpmcmc <y_file> <tobs_file> <Nt>
add_executable(run_ngpmcmc run_ngpmcmc.c)
target_link_libraries(run_ngpmcmc PRIVATE lafp_shared)

install(TARGETS run_ngpmcmc
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
