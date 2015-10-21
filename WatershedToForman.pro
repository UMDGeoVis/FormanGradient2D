TEMPLATE = app

CONFIG -= app_bundle
CONFIG -= qt

# Directories
DESTDIR = dist/
OBJECTS_DIR = build/

TARGET = watershedToForman

QMAKE_CXXFLAGS_RELEASE += -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++11

mac : {
QMAKE_CXXFLAGS += -stdlib=libc++
}

QMAKE_LFLAGS += -lstdc++

SOURCES += \
    source/main.cpp \
    source/LibMesh/Reader.cpp \
    source/LibMesh/Triangle.cpp \
    source/LibMesh/Vertex3D.cpp \
    source/LibForman/forman_feature.cpp \
    source/LibForman/formangradientvector.cpp \
    source/LibForman/io.cpp \
    source/LibForman/homotopy_expansion.cpp \
    source/LibMesh/mesh.cpp \
    source/LibForman/watershedtoforman.cpp \
    source/Performances/Timer.cpp \
    source/LibMesh/scalarfield.cpp

HEADERS += \
    source/LibMesh/Edge.h \
    source/LibMesh/Mesh.h \
    source/LibMesh/Reader.h \
    source/LibMesh/Sorting.h \
    source/LibMesh/Triangle.h \
    source/LibMesh/Vertex3D.h \
    source/LibForman/aux_forman_gradient_vector.h \
    source/LibForman/forman_arrow.h \
    source/LibForman/formangradientvector.h \
    source/WatershedAlgs/simulatedimmersion.h \
    source/Performances/Timer.h \
    source/Performances/Usage.h


