# Instruction to build libHaru from internal package source
INCLUDE(ExternalProject)
SET(EXTRA "")
IF(APPLE)
   SET(EXTRA  -DCMAKE_INSTALL_NAME_DIR=${CMAKE_INSTALL_PREFIX}/lib)
ENDIF(APPLE)

IF(NOT STANDALONE)
  SET(EXTRA ${EXTRA} -DBUILD_SHARED_LIBS=ON)
  SET(EXTRA ${EXTRA} -DBUILD_STATIC_LIBS=OFF)
ELSE(NOT STANDALONE)
	SET(EXTRA ${EXTRA} -DBUILD_STATIC_LIBS=ON)
	SET(EXTRA ${EXTRA} -DBUILD_SHARED_LIBS=OFF)
ENDIF(NOT STANDALONE) 

IF(CMAKE_TOOLCHAIN_FILE)
	SET(EXTRA ${EXTRA} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE})
ENDIF(CMAKE_TOOLCHAIN_FILE)

EXTERNALPROJECT_ADD(libHDF5
		#SOURCE_DIR ${PROJECT_SOURCE_DIR}/extern/hdf5-1.8.15-patch1
		URL https://www.hdfgroup.org/package/source-bzip2/?wpdmdl=4300&refresh=5a1c163c7f3491511790140
		PREFIX ${EXTERNAL_INSTALL_LOCATION}/libHDF5
		CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/libHDF5
			   -DBUILD_TESTING=OFF
			   -DBUILD_STATIC_EXECS=OFF
			   -DHDF5_BUILD_CPP_LIB=OFF
			   -DHDF5_BUILD_EXAMPLES=OFF
			   -DHDF5_BUILD_HL_LIB=OFF
			   -DHDF5_BUILD_TOOLS=OFF
			   -DHDF5_DISABLE_COMPILER_WARNINGS=ON
			   -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
			   -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
			   -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
			   ${EXTRA}
		INSTALL_COMMAND ""
		LOG_DOWNLOAD OFF
)
		
ExternalProject_Get_Property(libHDF5 install_dir)
include_directories(BEFORE SYSTEM "${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5/src;${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5-build/" )

MESSAGE( STATUS "libHDF5 will be built in ${install_dir}" )

#SET( hdf5_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR}/libHDF5/lib PATH)
# Set uncached variables as per standard.
# 	  set(hdf5_FOUND ON)
SET(HDF5_C_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5/src;${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5-build/)
MESSAGE(STATUS "HFF5 include dir is  ${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5/src")

IF(STANDALONE)
  SET(HDF5_C_LIBRARIES ${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5-build/bin/libhdf5-static${CMAKE_STATIC_LIBRARY_SUFFIX};-lz)
ELSE(STANDALONE)
  GET_FILENAME_COMPONENT(XYZ ${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5-build/bin/libhdf5-shared${CMAKE_SHARED_LIBRARY_SUFFIX} REALPATH)
  GET_FILENAME_COMPONENT( XYZ2 ${XYZ} REALPATH)
  SET(HDF5_C_LIBRARIES ${XYZ};-lz)
  
  IF(APPLE)
    SET(HDF5_PATTERN libhdf5-shared*${CMAKE_SHARED_LIBRARY_SUFFIX})
  ELSE(APPLE)
    SET(HDF5_PATTERN libhdf5-shared${CMAKE_SHARED_LIBRARY_SUFFIX}*)
  ENDIF(APPLE)

  INSTALL(DIRECTORY ${EXTERNAL_INSTALL_LOCATION}/libHDF5/src/libHDF5-build/bin/
          DESTINATION lib
          FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
	                   GROUP_EXECUTE GROUP_READ 
	                   WORLD_EXECUTE WORLD_READ
          FILES_MATCHING PATTERN ${HDF5_PATTERN}
 	)
ENDIF(STANDALONE)

MARK_AS_ADVANCED( HDF5_C_LIBRARIES HDF5_C_INCLUDE_DIR)
