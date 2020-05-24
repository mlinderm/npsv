set(SEQLIB_ROOT "${CMAKE_SOURCE_DIR}/lib/seqlib")

find_path(
    SEQLIB_INCLUDE_DIR
    NAMES "SeqLib/SeqLibCommon.h"
    PATHS ${SEQLIB_ROOT}
    PATH_SUFFIXES "src/seqlib"
)

find_library(
    SEQLIB_LIBRARY
    NAMES "libseqlib.a"
    PATHS ${SEQLIB_ROOT}
    PATH_SUFFIXES "src/seqlib/src"
)

find_library(
    HTS_LIBRARY
    NAMES "libhts.a"
    PATHS ${SEQLIB_ROOT}
    PATH_SUFFIXES "src/seqlib/htslib"
    NO_DEFAULT_PATH
)

find_library(
    BWA_LIBRARY
    NAMES "libbwa.a"
    PATHS ${SEQLIB_ROOT}
    PATH_SUFFIXES "src/seqlib/bwa"
    NO_DEFAULT_PATH
)


if ((EXISTS "${SEQLIB_INCLUDE_DIR}") AND (EXISTS "${SEQLIB_LIBRARY}"))
    set(SEQLIB_INCLUDE_DIRS ${SEQLIB_INCLUDE_DIR} ${SEQLIB_INCLUDE_DIR}/htslib)
    set(SEQLIB_LIBRARIES 
        ${SEQLIB_LIBRARY}
        ${HTS_LIBRARY}
        ${BWA_LIBRARY}
    )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(
        seqlib
        DEFAULT_MSG
        SEQLIB_INCLUDE_DIRS
        SEQLIB_LIBRARIES
    )
else()
    include(ExternalProject)
    # C++ language version needs to match main project
    ExternalProject_Add(
        seqlib
        SOURCE_DIR ${SEQLIB_ROOT}
        CONFIGURE_COMMAND ./configure
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make CXXFLAGS='-std=c++14'
        INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(seqlib source_dir binary_dir)
    set(SEQLIB_INCLUDE_DIRS ${source_dir} ${source_dir}/htslib)

    add_library(libseqlib IMPORTED STATIC GLOBAL)
    add_library(libhts IMPORTED STATIC GLOBAL)
    add_library(libbwa IMPORTED STATIC GLOBAL)

    add_dependencies(libseqlib seqlib)
    add_dependencies(libhts seqlib)
    add_dependencies(libbwa seqlib)

    set_target_properties(libseqlib PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/src/libseqlib.a"
    )

    set_target_properties(libhts PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/htslib/libhts.a"
    )

    set_target_properties(libbwa PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/bwa/libbwa.a"
    )

    set(SEQLIB_LIBRARIES libseqlib libhts libbwa)
endif()

if(EXISTS "${SEQLIB_INCLUDE_DIRS}")
    set(SEQLIB_FOUND 1)
else()
    set(SEQLIB_FOUND 0)
endif()