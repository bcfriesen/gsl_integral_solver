SET (gsl_integral_solver_SOURCES
    aux_g.c
    fred2_interp.c
    fred2_solver.c
    kernel.c
    main.c
    )

INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/src
    $ENV{GSL_DIR}/include
    )

LINK_DIRECTORIES(
    $ENV{GSL_DIR}/lib
    )

ADD_EXECUTABLE (${PROJECT_NAME} ${gsl_integral_solver_SOURCES})
TARGET_LINK_LIBRARIES (${PROJECT_NAME} gsl gslcblas)
