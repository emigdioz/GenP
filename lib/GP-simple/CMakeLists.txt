cmake_minimum_required(VERSION 3.0)

include(GNUInstallDirs)
add_library(GP-simple SHARED
    src/gpsingle.cpp
    src/paramssingle.cpp
    src/primitivessingle.cpp
    ../common/util/CmdLineParser.cc
)
include_directories(../common ../common/util)
set(HEADER_FILES
    src/gpsingle.h
    src/paramssingle.h
    src/primitivessingle.h
    ../common/Exception.h
    ../common/util/CmdLineException.h
    ../common/util/CmdLineParser.h
    ../common/util/Random.h
    ../common/util/Util.h)

set_target_properties(GP-simple PROPERTIES 
    VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER "${HEADER_FILES}")

configure_file(GP-simple.pc.in GP-simple.pc @ONLY)
target_include_directories(GP-simple PRIVATE .)
install(TARGETS GP-simple
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/lib/GP-simple/GP-simple.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)
