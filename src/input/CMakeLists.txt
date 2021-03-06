SET(INPUT_SOURCES "")
################################ FASTA ##################################
IF(PRF_INPUT_FASTA)
	SET(INPUT_SOURCES ${INPUT_SOURCES};FASTA/fasta.c)
ENDIF(PRF_INPUT_FASTA)

################################ PacBio #################################
IF(PRF_INPUT_HDF5)
	SET(INPUT_SOURCES ${INPUT_SOURCES};PacBio/pb_hdf5.c)
ENDIF(PRF_INPUT_HDF5)
IF(PRF_INPUT_PBBAM)
	SET(INPUT_SOURCES ${INPUT_SOURCES};PacBio/pb_bam.cpp)
	SET_PROPERTY(SOURCE PacBio/pb_bam.cpp APPEND PROPERTY COMPILE_FLAGS "-std=gnu++11 -I${Boost_DIR}")
ENDIF(PRF_INPUT_PBBAM)
IF(PRF_INPUT_HDF5 OR PRF_INPUT_PBBAM)
	SET(INPUT_SOURCES ${INPUT_SOURCES};PacBio/pb_common.c)
ENDIF(PRF_INPUT_HDF5 OR PRF_INPUT_PBBAM)

################################ Random #################################
# ADD_SUBDIRECTORY(Random)

#########################################################################
IF(NOT "${INPUT_SOURCES}" STREQUAL "")
	IF(BUILD_STATIC_LIBS)
		ADD_LIBRARY(INPUT_FORMAT_STATIC OBJECT ${INPUT_SOURCES})
	ENDIF(BUILD_STATIC_LIBS)
	IF(BUILD_SHARED_LIBS)
		ADD_LIBRARY(INPUT_FORMAT_SHARED OBJECT ${INPUT_SOURCES})
		SET_PROPERTY(TARGET INPUT_FORMAT_SHARED PROPERTY POSITION_INDEPENDENT_CODE ON)
	ENDIF(BUILD_SHARED_LIBS)

	IF(PRF_INPUT_PBBAM)
		IF(BUILD_STATIC_LIBS)
			SET_PROPERTY(TARGET INPUT_FORMAT_STATIC 
				     APPEND PROPERTY INCLUDE_DIRECTORIES "${PBBAM_INCLUDE_DIRS};${HTSLIB_DIR}")
			ADD_DEPENDENCIES(INPUT_FORMAT_STATIC PacBioBAM)
		ENDIF(BUILD_STATIC_LIBS)
		IF(BUILD_SHARED_LIBS)
			SET_PROPERTY(TARGET INPUT_FORMAT_SHARED
				     APPEND PROPERTY INCLUDE_DIRECTORIES "${PBBAM_INCLUDE_DIRS};${HTSLIB_DIR}")
			ADD_DEPENDENCIES(INPUT_FORMAT_SHARED PacBioBAM)
		ENDIF(BUILD_SHARED_LIBS)
	ENDIF(PRF_INPUT_PBBAM)

	IF(PRF_INPUT_HDF5)
		IF(BUILD_STATIC_LIBS)
			SET_PROPERTY(TARGET INPUT_FORMAT_STATIC 
									APPEND PROPERTY INCLUDE_DIRECTORIES  ${HDF5_C_INCLUDE_DIR})
			SET_PROPERTY(TARGET INPUT_FORMAT_STATIC 
									APPEND PROPERTY LINK_DEPENDS ${HDF5_C_LIBRARIES})
			IF(BUILD_HDF5) 
				ADD_DEPENDENCIES(INPUT_FORMAT_STATIC libHDF5)
			ENDIF(BUILD_HDF5) 
		ENDIF(BUILD_STATIC_LIBS)
		IF(BUILD_SHARED_LIBS)
			SET_PROPERTY(TARGET INPUT_FORMAT_SHARED
								APPEND PROPERTY INCLUDE_DIRECTORIES  ${HDF5_C_INCLUDE_DIR})
			SET_PROPERTY(TARGET INPUT_FORMAT_SHARED 
								APPEND PROPERTY LINK_DEPENDS ${HDF5_C_LIBRARIES})
			IF(BUILD_HDF5) 
				ADD_DEPENDENCIES(INPUT_FORMAT_SHARED libHDF5)
			ENDIF(BUILD_HDF5)
		ENDIF(BUILD_SHARED_LIBS)
	ENDIF(PRF_INPUT_HDF5)
ENDIF(NOT "${INPUT_SOURCES}" STREQUAL "")
