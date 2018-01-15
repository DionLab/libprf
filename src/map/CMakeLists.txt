IF(PRF_CORE_MAP)
	SET(JSON_C_SOURCES
		../../third-party/json-c/arraylist.c
		../../third-party/json-c/debug.c
		../../third-party/json-c/json_object.c
		../../third-party/json-c/json_pointer.c
		../../third-party/json-c/json_tokener.c
		../../third-party/json-c/json_util.c
		../../third-party/json-c/linkhash.c
		../../third-party/json-c/printbuf.c
		../../third-party/json-c/random_seed.c
	)
	SET_PROPERTY(SOURCE ${JSON_C_SOURCES} APPEND PROPERTY COMPILE_DEFINITIONS _GNU_SOURCE)
	ADD_LIBRARY(MAP OBJECT ${JSON_C_SOURCES} fixed_sse41.c map.c)
	SET_PROPERTY(SOURCE fixed_sse41.c APPEND PROPERTY COMPILE_FLAGS ${CMAKE_C_SSE41_FLAGS})
	TARGET_INCLUDE_DIRECTORIES(MAP PUBLIC ../../third-party/json-c)
	IF(BUILD_SHARED_LIBS)
		SET_PROPERTY(TARGET MAP PROPERTY POSITION_INDEPENDENT_CODE ON)
	ENDIF(BUILD_SHARED_LIBS)
ENDIF(PRF_CORE_MAP)
