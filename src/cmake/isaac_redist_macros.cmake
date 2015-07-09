################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## GNU GENERAL PUBLIC LICENSE Version 3
##
## You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
## along with this program. If not, see
## <https://github.com/illumina/licenses/>.
##
## file isaac_redist_macrocs.cmake 
##
## Configuration file for libxml2, libxslt library search, and redist
##
## author David Kimmel
##
################################################################################

##############################################################################
## Find a library with specific version, search on system defaults
macro(find_package_version libname version)
  string(TOUPPER ${libname} ${libname}_UPPER)

  find_package("${libname}" ${version})
 
  if(${${libname}_UPPER}_FOUND)
    if("${${${libname}_UPPER}_VERSION_STRING}" STREQUAL "${version}")
       message(" Found: ${libname}, correct version ${version}")
       message("   ${${libname}_UPPER}_INCLUDE_DIR = ${${${libname}_UPPER}_INCLUDE_DIR}")
       message("   ${${libname}_UPPER}_LIBRARIES = ${${${libname}_UPPER}_LIBRARIES}")
    else("${${${libname}_UPPER}_VERSION_STRING}" STREQUAL "${version}")
       message(" Not found: ${libname}, incorrect version ( ${${${libname}_UPPER}_VERSION} )")
       set(${${libname}_UPPER}_FOUND "FALSE")
    endif("${${${libname}_UPPER}_VERSION_STRING}" STREQUAL "${version}")
  endif(${${libname}_UPPER}_FOUND)

endmacro(find_package_version libname version)

##############################################################################
## Redist if not found (untar, configure, make install)
macro(redist_package name version args)

     if (NOT CMAKE_VERBOSE_MAKEFILE)
       list(APPEND REDIST_QUIET "OUTPUT_QUIET")
       list(APPEND REDIST_QUIET "ERROR_QUIET")
     endif (NOT CMAKE_VERBOSE_MAKEFILE)

     unset(HAVE_${name} CACHE)
     set(${name}_VERREQ "${name}-${version}")
     string(TOLOWER "${${name}_VERREQ}" nameb)
     set(tgzname "${${name}_REDIST_DIR}/${nameb}.tar.gz")
     message("    unpacking (${tgzname})") 
     execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${tgzname} RESULT_VARIABLE TMP_RESULT ${REDIST_QUIET})
     if (NOT TMP_RESULT)
        message(STATUS "Successfuly unpacked ${tgzname} from the distribution package...")
     else (NOT TMP_RESULT)
        message (FATAL_ERROR "Failed to unpack ${tgzname}")
     endif (NOT TMP_RESULT)
     
     message("    configuring ${name}")
     execute_process(COMMAND ./configure ${args} WORKING_DIRECTORY ${nameb} RESULT_VARIABLE TMP_RESULT ${REDIST_QUIET})
     if (NOT TMP_RESULT)
        message(STATUS "Successfuly configured ${nameb} from the distribution package...")
     else (NOT TMP_RESULT)
        message (FATAL_ERROR "Failed to configure ${nameb}")
     endif (NOT TMP_RESULT)
     
     message("    installing ${name}")
     execute_process(COMMAND make -j ${CMAKE_PARALLEL} install WORKING_DIRECTORY ${nameb} RESULT_VARIABLE TMP_RESULT ${REDIST_QUIET})
     if (NOT TMP_RESULT)
        message(STATUS "Successfuly built ${nameb} from the distribution package...")
     else (NOT TMP_RESULT)
        message (FATAL_ERROR "Failed to build ${name}")
     endif (NOT TMP_RESULT)
endmacro(redist_package name args)

##############################################################################
## Find a library given a path hint, assume version will be correct
macro(find_library_redist name pathhint header library)
    unset(${name}_LIBRARIES CACHE)
    # Search for library
    unset(${name}_LIBRARIES CACHE)
    find_library(${name}_LIBRARIES NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}${library}${CMAKE_STATIC_LIBRARY_SUFFIX} HINTS ${pathhint}/lib NO_DEFAULT_PATH)
    
    message(STATUS "Find library redist ${namenolib}, HINTS ${pathhint}/lib: ${${name}_LIBRARIES}")
    # Search for include path
    unset(${name}_INCLUDE_DIR CACHE)
    string(TOLOWER ${name} namel)
    find_path(${name}_INCLUDE_DIR ${header} HINTS ${pathhint}/include PATH_SUFFIXES ${namel} NO_DEFAULT_PATH)
    set(${name}_INCLUDE_DIR ${${name}_INCLUDE_DIR} CACHE STRING "lib BOOL" FORCE)

    if(${name}_INCLUDE_DIR AND ${name}_LIBRARIES)
        set (HAVE_${name} true CACHE BOOL "lib bool" FORCE)
        message (STATUS "Found redist ${name}  header: ${${name}_INCLUDE_DIR}/${header}")
        message (STATUS "Found redist ${name} library: ${${name}_LIBRARY}")
    endif(${name}_INCLUDE_DIR AND ${name}_LIBRARIES)

endmacro(find_library_redist name pathhint header)

