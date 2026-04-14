cmake_minimum_required(VERSION 3.14)
project(PionCloud VERSION 1.3 LANGUAGES CXX Fortran)

# ============================================================
#  Build options
# ============================================================
option(USE_ROOT   "Build with ROOT output support"    ON)
option(USE_CTEQ6  "Build with CTEQ6 Fortran PDF library" ON)
option(BUILD_TESTS "Build unit tests"                 ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release build by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_FLAGS_RELEASE     "-O3 -march=native -ffast-math")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")

# ============================================================
#  Sources
# ============================================================
set(LIB_SOURCES
    src/Calculator.cc
    src/OutputWriter.cc
    src/InputReader.cc
    src/pioncloud_c.cc
)

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ============================================================
#  Optional CTEQ6 Fortran library
#
#  Source: CTEQ6/Cteq6Pdf-2007.f  (v6.52, April 2007)
#  Tables: CTEQ6/cteq6m.tbl  (iset=1, NLO MSbar -- default)
#          CTEQ6/cteq6d.tbl  (iset=2, NLO DIS)
#          CTEQ6/cteq6l.tbl  (iset=3, LO)
#          CTEQ6/cteq6l1.tbl (iset=4, LO alphaS fitted)
#
#  CMake compiles Cteq6Pdf-2007.f as a Fortran translation unit
#  and links it into libPionCloudLib.a.  The C++ interface is
#  entirely in include/CTEQ6.hh via extern "C" bindings.
#
#  Table files are copied to the build directory at configure
#  time so the executable can find them at runtime.
#
#  gfortran name mangling: SETCTQ6 -> setctq6_, CTQPDF -> ctq6pdf_
#  This matches the extern "C" declarations in include/CTEQ6.hh.
# ============================================================
set(CTEQ6_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CTEQ6)
set(CTEQ6_SRC ${CTEQ6_DIR}/Cteq6Pdf-2007.f)

if(USE_CTEQ6)
    if(EXISTS ${CTEQ6_SRC})
        message(STATUS "CTEQ6 source found: ${CTEQ6_SRC}")
        list(APPEND LIB_SOURCES ${CTEQ6_SRC})
        set(CTEQ6_FOUND TRUE)
    else()
        message(WARNING
            "USE_CTEQ6=ON but ${CTEQ6_SRC} not found.\n"
            "Place Cteq6Pdf-2007.f in the CTEQ6/ directory and re-run cmake.\n"
            "Building without CTEQ6 -- F2n and ratio will be zero.")
        set(CTEQ6_FOUND FALSE)
    endif()
else()
    message(STATUS "CTEQ6 support disabled (USE_CTEQ6=OFF)")
    set(CTEQ6_FOUND FALSE)
endif()

# ============================================================
#  Library target
# ============================================================
add_library(PionCloudLib STATIC ${LIB_SOURCES})
target_include_directories(PionCloudLib PUBLIC ${INCLUDE_DIR})

if(CTEQ6_FOUND)
    target_compile_definitions(PionCloudLib PUBLIC HAVE_CTEQ6)
    # The CTEQ6 Fortran code reads its grid tables from files
    # at runtime using a relative path.  Copy (or symlink) the
    # table files into the build directory so the executable
    # finds them when run from there.
    file(GLOB CTEQ6_DATA_FILES "${CTEQ6_DIR}/*.tbl" "${CTEQ6_DIR}/*.dat")
    if(CTEQ6_DATA_FILES)
        foreach(_f ${CTEQ6_DATA_FILES})
            get_filename_component(_fname ${_f} NAME)
            configure_file(${_f}
                ${CMAKE_CURRENT_BINARY_DIR}/${_fname} COPYONLY)
        endforeach()
        message(STATUS "CTEQ6 data files copied to build directory.")
    else()
        message(WARNING
            "No CTEQ6 table files (*.tbl / *.dat) found in ${CTEQ6_DIR}.\n"
            "The executable will likely crash at runtime when CTEQ6::init()\n"
            "is called.  Place the grid files (e.g. cteq6l.tbl) in CTEQ6/.")
    endif()
endif()

# ============================================================
#  Optional ROOT support
#
#  Components used by OutputWriter.cc:
#    Core  -- TObject, TNamed, TROOT
#    RIO   -- TFile
#    Hist  -- TH1D
#    Graf  -- TGraph, TAxis, TLatex, TLegend
#    Gpad  -- TCanvas, TStyle
# ============================================================
if(USE_ROOT)
    find_package(ROOT QUIET
        COMPONENTS Core RIO Hist Graf Gpad)
    if(ROOT_FOUND)
        message(STATUS "ROOT found: ${ROOT_VERSION} at ${ROOT_INCLUDE_DIRS}")
        include(${ROOT_USE_FILE})
        target_compile_definitions(PionCloudLib PUBLIC HAVE_ROOT)
        target_include_directories(PionCloudLib PUBLIC ${ROOT_INCLUDE_DIRS})
        if(TARGET ROOT::Core)
            target_link_libraries(PionCloudLib PUBLIC
                ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad)
        else()
            target_link_libraries(PionCloudLib PUBLIC ${ROOT_LIBRARIES})
        endif()
    else()
        message(WARNING
            "ROOT not found -- building without ROOT output support.\n"
            "Install ROOT and re-run cmake, or disable with -DUSE_ROOT=OFF")
    endif()
else()
    message(STATUS "ROOT support disabled (USE_ROOT=OFF)")
endif()

# ============================================================
#  Executable
# ============================================================
add_executable(PionCloud src/main.cc)
target_link_libraries(PionCloud PRIVATE PionCloudLib)

# ============================================================
#  Optional tests
# ============================================================
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ============================================================
#  Install
# ============================================================
install(TARGETS PionCloud DESTINATION bin)
install(TARGETS PionCloudLib DESTINATION lib)
install(DIRECTORY include/ DESTINATION include/PionCloud)
