# Instruction to build PLplot from internal package source
SET(EXTRA "")
IF(APPLE)
   SET(EXTRA  -DCMAKE_INSTALL_NAME_DIR=${CMAKE_INSTALL_PREFIX}/lib)
ENDIF(APPLE)

INCLUDE(ExternalProject)
IF(STANDALONE)
	SET(SHARED_LIBS OFF)
ELSE(STANDALONE)
	SET(SHARED_LIBS ON)
ENDIF(STANDALONE)

IF(NOT BUILD_GD)
	EXTERNALPROJECT_ADD(PLplot
			URL ${PROJECT_SOURCE_DIR}/third-party/plplot-5.11.1.tar.gz
			PREFIX ${EXTERNAL_INSTALL_LOCATION}/PLplot
			CMAKE_ARGS 	-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/PLplot
			-DBUILD_DOC=OFF
			-DBUILD_TEST=OFF
			-DBUILD_DOX_DOC=OFF
			-DBUILD_SHARED_LIBS=${SHARED_LIBS}
			-DDEFAULT_NO_CAIRO_DEVICES=ON
			-DDEFAULT_NO_QT_DEVICES=ON
			-DDEFAULT_NO_BINDINGS=ON
			-DDEFAULT_NO_DEVICES=ON
			-DENABLE_DYNDRIVERS=OFF
			-DPLD_jpeg=ON
			-DPLD_png=ON
			${EXTRA}
	)
ELSE(NOT BUILD_GD)
	EXTERNALPROJECT_ADD(PLplot
			DEPENDS libGD
			URL ${PROJECT_SOURCE_DIR}/third-party/plplot-5.11.1.tar.gz
			PREFIX ${EXTERNAL_INSTALL_LOCATION}/PLplot
			CMAKE_ARGS 	-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/PLplot
			-DBUILD_DOC=OFF
			-DBUILD_TEST=OFF
			-DBUILD_DOX_DOC=OFF
			-DBUILD_SHARED_LIBS=${SHARED_LIBS}
			-DDEFAULT_NO_CAIRO_DEVICES=ON
			-DDEFAULT_NO_QT_DEVICES=ON
			-DDEFAULT_NO_BINDINGS=ON
			-DDEFAULT_NO_DEVICES=ON
			-DENABLE_DYNDRIVERS=OFF
			-DPLD_jpeg=ON
			-DPLD_png=ON
			-DGD_INCLUDE_DIR:PATH=${GD_INCLUDE_DIRS}
			-DGD_LIBRARY:PATH=${GD_LIBRARIES}
			-DUSE_RPATH=ON
			${EXTRA}
	)
ENDIF(NOT BUILD_GD)

ExternalProject_Get_Property(PLplot install_dir)
#include_directories(BEFORE SYSTEM  ${EXTERNAL_INSTALL_LOCATION}/PLplot/include)
MESSAGE(STATUS "PLplot will be built in ${install_dir}")

SET( PLplot_LIBRARY_DIRS ${install_dir}/src/PLplot-build/lib PATH)
# Set uncached variables as per standard.
# 	  set(hpdf_FOUND ON)
SET(PLplot_INCLUDE_DIRS ${install_dir}/src/PLplot/include)
MESSAGE(STATUS "PLplot include dir: ${PLplot_INCLUDE_DIRS}")

IF(STANDALONE)
	SET(PLplot_LIBRARIES 
				${install_dir}/src/PLplot-build/src/libplplot${CMAKE_STATIC_LIBRARY_SUFFIX};
				${install_dir}/src/PLplot-build/lib/csa/libcsirocsa${CMAKE_STATIC_LIBRARY_SUFFIX};
				${install_dir}/src/PLplot-build/lib/qsastime/libqsastime${CMAKE_STATIC_LIBRARY_SUFFIX};
				${GD_LIBRARIES}
	)
ELSE(STANDALONE)
	IF(NOT BUILD_GD)
		SET(PLplot_LIBRARIES 
					${install_dir}/src/PLplot-build/src/libplplot${CMAKE_SHARED_LIBRARY_SUFFIX};
					${install_dir}/src/PLplot-build/lib/csa/libcsirocsa${CMAKE_SHARED_LIBRARY_SUFFIX};
					${install_dir}/src/PLplot-build/lib/qsastime/libqsastime${CMAKE_SHARED_LIBRARY_SUFFIX}
		)
	ELSE(NOT BUILD_GD)
		SET(PLplot_LIBRARIES 
					${install_dir}/src/PLplot-build/src/libplplot${CMAKE_SHARED_LIBRARY_SUFFIX};
					${install_dir}/src/PLplot-build/lib/csa/libcsirocsa${CMAKE_SHARED_LIBRARY_SUFFIX};
					${install_dir}/src/PLplot-build/lib/qsastime/libqsastime${CMAKE_SHARED_LIBRARY_SUFFIX};
					${install_dir}/libGD/lib/libgd${CMAKE_SHARED_LIBRARY_SUFFIX}
		)
	ENDIF(NOT BUILD_GD)
	
	find_library(plplot1 NAMES plplot
	             PATHS "${install_dir}/src/PLplot-build/src/"
	             NO_DEFAULT_PATH )
	MESSAGE(STATUS "plplot is ${plplot1}")
	
ENDIF(STANDALONE)

MARK_AS_ADVANCED( PLplot_INCLUDE_DIRS PLplot_LIBRARIES)
