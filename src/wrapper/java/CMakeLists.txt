find_package(Java COMPONENTS Development)
find_package(JNI REQUIRED)
include(UseJava)
set(CMAKE_JAVA_TARGET_OUTPUT_DIR ${CMAKE_BINARY_DIR}/java)
MESSAGE(STATUS "Java build in ${CMAKE_JAVA_TARGET_OUTPUT_DIR}")
# set(CMAKE_JAVA_COMPILE_FLAGS "-source" "1.6" "-target" "1.6")
# 
include_directories(${JNI_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_JNI_TARGET TRUE)
set(java_sources Profile.java Matrix.java)
add_jar(PfTools ${java_sources})
Create_javah(
	TARGET jni_headers
# 	GENERATED_FILES jni_headers
  CLASSES PfTools.Profile PfTools.Matrix
  CLASSPATH PfTools
# 	CLASSPATH ${CMAKE_CURRENT_SOURCE_DIR}
   OUTPUT_DIR ${CMAKE_JAVA_TARGET_OUTPUT_DIR}
)
add_dependencies(jni_headers PfTools)
add_library(JNI OBJECT jniwrapper.c)
set_property(TARGET JNI PROPERTY POSITION_INDEPENDENT_CODE ON)
add_dependencies(JNI jni_headers)

# MESSAGE(STATUS "JNI header ${jni_headers}")

# install_jar(PfTools ${CMAKE_JAVA_TARGET_OUTPUT_DIR})
# install_jni_symlink(prf_shared .)
# 
# get_target_property(_jarFile HelloWorld JAR_FILE)
# get_target_property(_classDir HelloWorld CLASSDIR)
# 
# message(STATUS "Jar file ${_jarFile}")
# message(STATUS "Class compiled to ${_classDir}")
# 
# add_test(NAME TestHelloWorld COMMAND ${Java_JAVA_EXECUTABLE} -cp ${_jarFile} HelloWorld)
