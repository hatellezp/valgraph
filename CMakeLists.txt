cmake_minimum_required(VERSION 3.16)
project(valgraphcore VERSION 1.0)


### add subdirectories ###
add_subdirectory(core)

### adding executables and lib ###
add_executable(main main.c)
# add_executable(poly_coeff poly_coeff.c)
# add_executable(eval_time eval_time.c vgraph.c)
# add_executable(test_polynomial test_polynomial.c)

add_library(vgraph STATIC vgraph.c)


### include directories ###
# include_directories(core)


### target link libraries ###
target_link_libraries(vgraph PUBLIC core fftw3 m gsl)

# target_link_libraries(poly_coeff PUBLIC core fftw3 m gsl)
# target_link_libraries(eval_time PUBLIC core fftw3 m gsl)
target_link_libraries(main PUBLIC core gsl m)
# target_link_libraries(test_polynomial PUBLIC core gsl m)


# target include directories
target_include_directories(
        main PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/core"
)


set(CMAKE_C_STANDARD 90)
