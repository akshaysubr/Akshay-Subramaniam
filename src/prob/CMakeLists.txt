set(DRIVER_PATH ../driver)
add_executable(pp_allprocs pp_allprocs.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(pp_lowmem pp_lowmem.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
