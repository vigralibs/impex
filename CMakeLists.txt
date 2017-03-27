cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

enable_testing()

project(vigra_imageio VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/vigra_cmake")
include(VigraAddDep)

vigra_add_dep(vigra_core LIVE)

option(ENABLE_PNG "Build with PNG support.")
option(ENABLE_OPENEXR "Build with OpenEXR support.")
option(ENABLE_JPEG "Build with JPEG support.")
option(ENABLE_TIFF "Build with TIFF support.")
option(VIGRA_IMAGEIO_STATIC_LIB "Build a static lib." OFF)

if(ENABLE_PNG)
  vigra_add_dep(PNG)
endif()

if(ENABLE_OPENEXR)
  vigra_add_dep(OpenEXR)
endif()

if(ENABLE_JPEG)
  vigra_add_dep(JPEG)
endif()

if(ENABLE_TIFF)
  vigra_add_dep(TIFF)
endif()

IF(VIGRA_IMAGEIO_STATIC_LIB)
    SET(LIBTYPE STATIC)
ELSE()
    SET(LIBTYPE SHARED)
ENDIF()

add_library(vigra_imageio ${LIBTYPE}
    src/bmp.cxx
    src/byteorder.cxx
    src/codecmanager.cxx
    src/compression.cxx
    src/exr.cxx
    src/gif.cxx
    src/hdr.cxx
    #src/hdf5impex.cxx
    #src/hdf5_rf_impex.cxx
    src/iccjpeg.c
    src/imageinfo.cxx
    src/jpeg.cxx
    src/lz4.c
    src/png.cxx
    src/pnm.cxx
    src/rgbe.c
    src/sifImport.cxx
    src/sun.cxx
    src/tiff.cxx
    src/viff.cxx
    src/void_vector.cxx
)

target_include_directories(vigra_imageio PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

target_link_libraries(vigra_imageio PUBLIC vigra_core)

if(MSVC)
  if(VIGRA_IMAGEIO_STATIC_LIB)
    target_compile_definitions(vigra_imageio PUBLIC VIGRA_STATIC_LIB)
  else()
    target_compile_definitions(vigra_imageio PRIVATE VIGRA_DLL)
  endif()
endif()

# set(SOVERSION 11)  # increment this after changing the vigraimpex library
# IF(MACOSX)
    # SET_TARGET_PROPERTIES(vigraimpex PROPERTIES VERSION ${SOVERSION}.${vigra_version}
                         # SOVERSION ${SOVERSION} INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
# ELSE()
    # SET_TARGET_PROPERTIES(vigraimpex PROPERTIES VERSION ${SOVERSION}.${vigra_version} SOVERSION ${SOVERSION})
# ENDIF()
MESSAGE(STATUS "FIXME: we must add the SO version if we build a shared library")

#target_link_libraries(vigra_imageio PUBLIC ZLIB::ZLIB)

if(ENABLE_PNG)
  target_compile_definitions(vigra_imageio PRIVATE HasPNG)
  target_link_libraries(vigra_imageio PUBLIC PNG::PNG)
endif()

if(ENABLE_OPENEXR)
  target_compile_definitions(vigra_imageio PRIVATE HasEXR)
  target_link_libraries(vigra_imageio PUBLIC OpenEXR::OpenEXR)
endif()

if(ENABLE_JPEG)
  target_compile_definitions(vigra_imageio PRIVATE HasJPEG)
  target_link_libraries(vigra_imageio PUBLIC JPEG::JPEG)
endif()

if(ENABLE_TIFF)
  target_compile_definitions(vigra_imageio PRIVATE HasTIFF)
  target_link_libraries(vigra_imageio PUBLIC TIFF::TIFF)
endif()

install(DIRECTORY include/ DESTINATION include)

add_subdirectory(test)
