set(DRIVER_PATH ../driver)
add_executable(postprocess postprocess.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/interfaces.f90
	${DRIVER_PATH}/functions.f90
)
