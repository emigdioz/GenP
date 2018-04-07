cmake_minimum_required(VERSION 3.0)

include(GNUInstallDirs)
add_library(GP-parallel SHARED
    src/gpparallel.cpp
    src/paramsparallel.cpp
    src/primitivesparallel.cpp
    src/CPU.cpp
    src/GPU.cpp
    ../common/util/CmdLineParser.cc
)
include_directories(../common ../common/util)
set(HEADER_FILES
    src/gpparallel.h
    src/paramsparallel.h
    src/primitivesparallel.h
    src/CPU.h
    src/GPU.h
    ../common/Exception.h
    ../common/util/CmdLineException.h
    ../common/util/CmdLineParser.h
    ../common/util/Random.h
    ../common/util/Util.h)

set(CL_FILES
    src/CL/cl_d3d10.h
    src/CL/cl_ext.h
    src/CL/cl_gl_ext.h
    src/CL/cl_gl.h
    src/CL/cl_platform.h
    src/CL/cl.h
    src/CL/cl.hpp
    src/CL/opencl.hpp)

set_target_properties(GP-parallel PROPERTIES 
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER "${HEADER_FILES}")
configure_file(GP-parallel.pc.in GP-parallel.pc @ONLY)
target_include_directories(GP-parallel PRIVATE CL)
target_link_libraries(GP-parallel -lOpenCL -lopenblas)
install(TARGETS GP-parallel
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/lib/GP-parallel/GP-parallel.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
