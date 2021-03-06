INCLUDE_DIRECTORIES(BEFORE SYSTEM "${PFTOOLS_SOURCE_DIR}/src/include/")

################################ PROFILES ####################################
SET(PROFILE_SOURCES
    io.c
    profile.c
    Normalization.c
)
IF(BUILD_STATIC_LIBS)
	ADD_LIBRARY(PROFILE_STATIC OBJECT ${PROFILE_SOURCES})
ENDIF(BUILD_STATIC_LIBS)
IF(BUILD_SHARED_LIBS)
	ADD_LIBRARY(PROFILE_SHARED OBJECT ${PROFILE_SOURCES})
	SET_PROPERTY(TARGET PROFILE_SHARED PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF(BUILD_SHARED_LIBS)
########################### REGULAR EXPRESSION ###############################
IF(PRF_CORE_PCRE)
	SET(REGEX_SOURCE pfregexp.c)
		ADD_LIBRARY(REGEX_STATIC OBJECT ${REGEX_SOURCE})
	SET_PROPERTY(SOURCE pfregexp.c APPEND PROPERTY INCLUDE_DIRECTORIES "${PCRE_INCLUDE_DIR}")
  IF(BUILD_SHARED_LIBS)
		ADD_LIBRARY(REGEX_SHARED OBJECT ${REGEX_SOURCE})
		SET_PROPERTY(TARGET REGEX_SHARED PROPERTY POSITION_INDEPENDENT_CODE ON)
	ENDIF(BUILD_SHARED_LIBS)
ENDIF(PRF_CORE_PCRE)

############################### STATISTICS ###################################
SET(STATISTICS_SOURCES
    statistics.c
)
SET_PROPERTY(SOURCE statistics.c APPEND PROPERTY COMPILE_FLAGS ${CMAKE_C_SSE41_FLAGS})
IF(BUILD_STATIC_LIBS)
	ADD_LIBRARY(STATISTICS_STATIC OBJECT ${STATISTICS_SOURCES})
ENDIF(BUILD_STATIC_LIBS)
IF(BUILD_SHARED_LIBS)
	ADD_LIBRARY(STATISTICS_SHARED OBJECT ${STATISTICS_SOURCES})
	SET_PROPERTY(TARGET STATISTICS_SHARED PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF(BUILD_SHARED_LIBS)
