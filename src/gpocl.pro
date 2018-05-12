TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lOpenCL
#LIBS += -static

SOURCES += \
    CPU.cc \
    GP.cc \
    gpocl.cc \
    GPU.cc \
    Params.cc \
    Primitives.cc \
    common/util/CmdLineParser.cc

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    CPU.h \
    GP.h \
    GPU.h \
    Params.h \
    Primitives.h \
    CL/cl.h \
    CL/cl.hpp \
    CL/cl_d3d10.h \
    CL/cl_ext.h \
    CL/cl_gl.h \
    CL/cl_gl_ext.h \
    CL/cl_platform.h \
    CL/opencl.h \
    common/util/CmdLineException.h \
    common/util/CmdLineParser.h \
    common/util/Random.h \
    common/util/Util.h \
    common/Exception.h

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
