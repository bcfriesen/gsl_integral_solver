SET (gsl_integral_solver_SOURCES
    main.c
    )

INCLUDE_DIRECTORIES(
    $ENV{GSL_DIR}/include
    )

LINK_DIRECTORIES(
    $ENV{GSL_DIR}/lib
    )

ADD_EXECUTABLE (${PROJECT_NAME} ${gsl_integral_solver_SOURCES})
TARGET_LINK_LIBRARIES (${PROJECT_NAME} gsl gslcblas)
