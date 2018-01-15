FIND_PROGRAM(PATCH_EXECUTABLE patch)
find_package(PythonInterp REQUIRED)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)
ExternalProject_Add(
    gtest
    PREFIX ${EXTERNAL_INSTALL_LOCATION}/gtest
#    GIT_REPOSITORY https://gitlab.isb-sib.ch/tschuepb/googletest.git
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG "release-1.8.0"
    TIMEOUT 10
#    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    PATCH_COMMAND ${PATCH_EXECUTABLE} -t -Np1 < ${CMAKE_CURRENT_SOURCE_DIR}/Patches/gtest.patch || true
    LOG_DOWNLOAD ON
)
ExternalProject_Get_Property(gtest source_dir)
MESSAGE(STATUS "Gtest is in ${source_dir}")
SET(GTEST_DIR ${source_dir})
function(pumpit inputfilename outputfilename)
	STRING(REPLACE ".pump" "" name "${inputfilename}")
	add_custom_command(OUTPUT "${outputfilename}"
		COMMAND python2.7 ${GTEST_DIR}/googletest/scripts/pump.py ${inputfilename}
		COMMAND ${CMAKE_COMMAND} -E rename "${name}" "${outputfilename}"
		DEPENDS "${inputfilename}" gtest
		COMMENT "Pumping ${inputfilename}"
	)
endfunction(pumpit)
