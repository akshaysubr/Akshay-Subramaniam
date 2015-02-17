set(DRIVER_PATH ../driver)
add_executable(vortstats_allprocs vortstats_allprocs.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(vortstats_lowmem vortstats_lowmem.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(vortz_YZint vortz_YZint.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(wall_torque wall_torque.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
add_executable(vort_eqterms vort_eqterms.f90
	${DRIVER_PATH}/IO.f90
	${DRIVER_PATH}/globals.f90
	${DRIVER_PATH}/operators.f90
	${DRIVER_PATH}/functions.f90
	${DRIVER_PATH}/interfaces.f90
)
