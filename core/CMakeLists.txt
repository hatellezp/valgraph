add_library(core matrix.c algrank.c time_eval.c gnuplot_i.c  csv_reader.c polynomial.c  polynomial_matrix.c  polynomial_vector.c todo.c c_error.c c_strerror.c c_stream_err.c polynomial_linalg.c lib.c test_lib.c test_lib.h)

install(TARGETS core DESTINATION lib)
install(FILES matrix.h algrank.h c_errno.h lib.h  DESTINATION include)