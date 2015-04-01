set(DRIVER_PATH ../driver)
add_executable(test_parallel test_parallel.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/mpi.f90
    ${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(vort_eqterms vort_eqterms.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/mpi.f90
    ${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(primitive_budgets primitive_budgets.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/mpi.f90
    ${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
