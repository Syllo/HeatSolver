add_executable(heatSolverResilienceTest main.c heatsolver.c gauss-seidel.c jacobi.c over-relaxation.c initialize.c)
target_include_directories(heatSolverResilienceTest PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(heatSolverResilienceTest PRIVATE m)
set_property(TARGET heatSolverResilienceTest
             PROPERTY C_STANDARD 11)

find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(heatSolverResilienceTest PRIVATE OpenMP::OpenMP_C)
endif()

# Compile Options
include(compile-flags-helpers)
include(${PROJECT_SOURCE_DIR}/optimization_flags.cmake)

if (DEFINED ADDITIONAL_BENCHMARK_COMPILE_OPTIONS)
  add_compiler_option_to_target_type(heatSolverResilienceTest Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_COMPILE_OPTIONS})
endif()

foreach(compile_type IN ITEMS Release RelWithDebInfo)
  add_compiler_option_to_target_type(heatSolverResilienceTest ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_COMPILE_OPTIONS})
  add_linker_option_to_target_type(heatSolverResilienceTest ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_LINK_OPTIONS})
endforeach()

add_compiler_option_to_target_type(heatSolverResilienceTest Debug PRIVATE ${ADDITIONAL_DEBUG_COMPILE_OPTIONS})

# Linker Options

if (DEFINED ADDITIONAL_BENCHMARK_LINK_OPTIONS)
  add_linker_option_to_target_type(heatSolverResilienceTest Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_LINK_OPTIONS})
endif()

add_sanitizers_to_target(heatSolverResilienceTest Debug PRIVATE address undefined)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if((result) AND USE_IPO)
  set_property(TARGET heatSolverResilienceTest PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

install(TARGETS heatSolverResilienceTest RUNTIME DESTINATION bin)
