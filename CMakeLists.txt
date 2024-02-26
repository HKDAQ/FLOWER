# Minimum cmake verison 3.1 required for the variable CMAKE_CXX_STANDARD
cmake_minimum_required(VERSION 3.1)


set(PROJECT_NAME FLOWER)

list(APPEND CMAKE_MODULE_PATH
        $ENV{HK_PILOT_DIR}/cmake
        )
include(GetVersionGit)
MESSAGE("Preparing project ${PROJECT_NAME} version ${PROJECT_VERSION}" )
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION})
include(HKPackageBuilder)
hkbuilder_prepare_project()
set( PUBLIC_EXT_LIBS )
hk_check_dependencies()

include_directories(BEFORE
        ${PROJECT_SOURCE_DIR}
        )

add_compile_definitions(HEMI_CUDA_DISABLE)
add_compile_definitions(NOSKLIBRARIES)

######################
#   src             #
######################

#add_library (fiTQunFortranLib sortzv.F)


# If there were other libraries in this package on which this library depends, then they would be put in this variable

# If there were other libraries in another package (ROOT,  on which this library depends, then they would be put in this variable
list ( APPEND PUBLIC_EXT_LIBS
        ROOT::Core
        ROOT::MathMore
        ROOT::RIO
        ROOT::RooFit
        WCSim::WCSimRoot)

set(sources
    WCSimFLOWER.cpp
        )


#set(exe_sources
#        runfiTQun.cc
#        )


set(headers
    WCSimFLOWER.h
        )

# ###########
# # Library #
# ###########


## WCSimFLOWERDict
#set (dict_headers
#        WCSimFLOWER.h
#        )
#set( linkdef_header
#        ${PROJECT_SOURCE_DIR}/WCSimFLOWER.h
#        )
#set(pcm_file ${CMAKE_CURRENT_BINARY_DIR}/libWCSimFLOWER_rdict.pcm)
#set(rootmap_file ${CMAKE_CURRENT_BINARY_DIR}/libWCSimFLOWER.rootmap)
#ROOT_GENERATE_DICTIONARY(G__WCSimFlower ${dict_headers} LINKDEF ${linkdef_header})
#pbuilder_install_files(${LIB_INSTALL_DIR} ${pcm_file} ${rootmap_file})
#
#set(sources ${sources} ${CMAKE_CURRENT_BINARY_DIR}/G__WCSimFLOWER.cxx)


pbuilder_library(
        TARGET FLOWER
        SOURCES ${sources}
        PROJECT_LIBRARIES ${PACKAGE_LIBS}
        PUBLIC_EXTERNAL_LIBRARIES ${PUBLIC_EXT_LIBS}
        PRIVATE_EXTERNAL_LIBRARIES ${PRIVATE_EXT_LIBS}
)

pbuilder_install_headers(${headers})


#########
# executable
#########
#
#pbuilder_executables(
#        SOURCES ${exe_sources}
#        TARGETS_VAR programs
#        PROJECT_LIBRARIES fiTQunLib fiTQunFortranLib
#        PUBLIC_EXTERNAL_LIBRARIES ${PUBLIC_EXT_LIBS}
#        PRIVATE_EXTERNAL_LIBRARIES ${PRIVATE_EXT_LIBS}
#)

#########
# target exportation
#########

pbuilder_component_install_and_export(
        COMPONENT Library #Executable
        LIBTARGETS FLOWER
#        EXETARGETS ${programs}
)