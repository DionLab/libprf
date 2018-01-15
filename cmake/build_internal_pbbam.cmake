# htslib
#FIND_PACKAGE(htslib QUIET)
IF (NOT htslib_FOUND)
	MESSAGE(STATUS "HTS lib is required to compile PBBAM, trying downloading")
	find_package(ZLIB REQUIRED)
	IF(APPLE)
		SET(EXTRA ${CMAKE_INSTALL_NAME_TOOL} -id @rpath/libhts.dylib ${EXTERNAL_INSTALL_LOCATION}/htslib/src/htslib/libhts.dylib)
		ExternalProject_Add(
				htslib
				PREFIX  ${EXTERNAL_INSTALL_LOCATION}/htslib
				GIT_REPOSITORY https://github.com/samtools/htslib.git
				TIMEOUT 10
				UPDATE_COMMAND ${GIT_EXECUTABLE} pull
				CONFIGURE_COMMAND autoreconf && CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ./configure --disable-lzma --disable-bz2 --disable-libcurl --disable-s3 --disable-gcs
				BUILD_COMMAND make -j 4
				BUILD_IN_SOURCE 1
				INSTALL_COMMAND ${EXTRA}
				LOG_DOWNLOAD OFF
		)
	ELSE(APPLE)
		ExternalProject_Add(
				htslib
				PREFIX  ${EXTERNAL_INSTALL_LOCATION}/htslib
				GIT_REPOSITORY https://github.com/samtools/htslib.git
				TIMEOUT 10
				UPDATE_COMMAND ${GIT_EXECUTABLE} pull
				CONFIGURE_COMMAND autoreconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl --disable-s3 --disable-gcs
				BUILD_COMMAND make -j 4
				BUILD_IN_SOURCE 1
				INSTALL_COMMAND ""
				LOG_DOWNLOAD OFF
		)
	ENDIF(APPLE)

	
	ExternalProject_Get_Property(htslib source_dir)
	MESSAGE(STATUS "htslib downloaded into ${source_dir}")
	SET(HTSLIB_DIR "${source_dir}")
	SET(HTSLIB_LIBRARIES "${source_dir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX};-lz")
	
  INSTALL(FILES "${source_dir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}"
          DESTINATION lib
          PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
											GROUP_EXECUTE GROUP_READ 
											WORLD_EXECUTE WORLD_READ
 	)
	
ENDIF(NOT htslib_FOUND)

# Boost version 1.6x or above 
FIND_PACKAGE(Boost )
IF ((NOT Boost_FOUND) OR ( ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} VERSION_LESS 1.60))
	MESSAGE(STATUS "Boost 1.6 or higher is required to compile PBBAM, trying downloading")
	ExternalProject_Add(
			Boost
			PREFIX  ${EXTERNAL_INSTALL_LOCATION}/Boost
			URL https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2
			TIMEOUT 60
			UPDATE_COMMAND ""
			CONFIGURE_COMMAND ""
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
			LOG_DOWNLOAD OFF
	)
	ExternalProject_Get_Property(Boost source_dir)
	MESSAGE(STATUS "Boost downloaded into ${source_dir}")
	SET(Boost_DIR "${source_dir}" CACHE PATH "Boost directory" FORCE)
	SET(Boost_INCLUDE_DIRS "${source_dir}" CACHE PATH "Boost directory" FORCE)
ENDIF ((NOT Boost_FOUND) OR ( ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}  VERSION_LESS 1.60))

IF(APPLE)
        SET(PB_RPATH -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=FALSE -DCMAKE_INSTALL_NAME_DIR=@rpath)
ELSE(APPLE)
        SET(PB_RPATH -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=FALSE -DCMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH})
ENDIF(APPLE)

IF(CMAKE_TOOLCHAIN_FILE)
	SET(EXTRA ${EXTRA} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE})
ELSE(CMAKE_TOOLCHAIN_FILE)
	SET(EXTRA "")
ENDIF(CMAKE_TOOLCHAIN_FILE)

ExternalProject_Add(
			PacBioBAM
			PREFIX ${EXTERNAL_INSTALL_LOCATION}/pbbam
			URL https://github.com/PacificBiosciences/pbbam/archive/0.19.0.tar.gz
# 			URL https://github.com/PacificBiosciences/pbbam/archive/0.23.0.tar.gz
			#GIT_REPOSITORY https://github.com/PacificBiosciences/pbbam.git
			#GIT_REPOSITORY https://gitlab.isb-sib.ch/tschuepb/pbbam.git
			#GIT_TAG "49619a0423c467f69c4ae423b04764eb9a5a67f5"
			TIMEOUT 60
			UPDATE_COMMAND ${GIT_EXECUTABLE} pull
			PATCH_COMMAND ${CMAKE_COMMAND} -E echo "INSTALL(TARGETS pbbam LIBRARY DESTINATION ${PROJECT_BINARY_DIR}/PacBioBAM)" >> ${EXTERNAL_INSTALL_LOCATION}/pbbam/src/PacBioBAM/src/CMakeLists.txt
			LOG_DOWNLOAD OFF
			CMAKE_ARGS -DPacBioBAM_build_tools=OFF
                        -DPacBioBAM_build_tests=OFF
                        -DPacBioBAM_build_docs=OFF
                        -DPacBioBAM_build_shared=ON
                        -DBoost_INCLUDE_DIRS:FILEPATH=${Boost_DIR}
                        -DHTSLIB_INCLUDE_DIRS=${HTSLIB_DIR}
                        -DHTSLIB_LIBRARIES=${HTSLIB_LIBRARIES}
                        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/PacBioBAM
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                        ${PB_RPATH} ${EXTRA}
)

ExternalProject_Get_Property(PacBioBAM install_dir)
GET_FILENAME_COMPONENT(XYZ ${PROJECT_BINARY_DIR}/PacBioBAM/libpbbam${CMAKE_SHARED_LIBRARY_SUFFIX} REALPATH)
INSTALL(FILES ${XYZ}
				DESTINATION lib
				PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
										GROUP_EXECUTE GROUP_READ 
										WORLD_EXECUTE WORLD_READ
)

IF ((NOT Boost_FOUND) OR ( ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}  VERSION_LESS 1.60))
	add_dependencies(PacBioBAM Boost)
	SET(PBBAM_INCLUDE_DIRS "${Boost_DIR};${install_dir}/src/PacBioBAM/include")
ELSE ((NOT Boost_FOUND) OR ( ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}  VERSION_LESS 1.60))
	SET(PBBAM_INCLUDE_DIRS "${install_dir}/src/PacBioBAM/include")
ENDIF ((NOT Boost_FOUND) OR ( ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}  VERSION_LESS 1.60))
IF (NOT htslib_FOUND)
	add_dependencies(PacBioBAM htslib)
ENDIF(NOT htslib_FOUND)

MESSAGE(STATUS "PBBam not found on the system, downloading from GitHub...")
MESSAGE(STATUS "PBBam built in ${install_dir}")


SET(PBBAM_LIBRARIES ${PROJECT_BINARY_DIR}/PacBioBAM/libpbbam${CMAKE_SHARED_LIBRARY_SUFFIX})
IF(NOT htslib_FOUND) 
	SET(PBBAM_LIBRARIES "${PBBAM_LIBRARIES};${HTSLIB_LIBRARIES}")
ENDIF(NOT htslib_FOUND) 
SET(PBBAM_LIBRARY_DIRS  ${PROJECT_BINARY_DIR}/PacBioBAM)
