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
################################################################################
##
## file macros.cmake
##
## CMake configuration file for common installation macros
##
## authors: Roman Petrovski, Mauricio Varea
##
################################################################################

macro(configure_files srcDir destDir pattern)
    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        message(STATUS "Configuring file ${srcDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} @ONLY IMMEDIATE)
    endforeach(templateFile)
endmacro(configure_files)

macro(install_files srcDir destDir pattern perm)
    file(GLOB templateFiles ${srcDir}/${pattern})
    file(INSTALL DESTINATION ${destDir} TYPE FILE
         FILES ${templateFiles} PERMISSIONS ${perm})
endmacro(install_files)

macro(configure_files_recursively srcDir destDir pattern)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        message(STATUS "Configuring file ${srcDir}/${templateFile}")
        configure_file(${srcDir}/${templateFile} ${destDir}/${templateFile} @ONLY IMMEDIATE)
    endforeach(templateFile)
endmacro(configure_files_recursively)

macro(install_files_recursively srcDir destDir pattern perm)
    file(GLOB_RECURSE templateFiles RELATIVE ${srcDir} ${srcDir}/${pattern})
    foreach(templateFile ${templateFiles})
        get_filename_component(DIRNAME "${templateFile}" PATH)
        file(INSTALL DESTINATION ${destDir}/${DIRNAME} TYPE FILE
             FILES ${srcDir}/${templateFile} PERMISSIONS ${perm} PATTERN CMakeLists.txt EXCLUDE )
    endforeach(templateFile)
endmacro(install_files_recursively)

#   
# Macro to find libraries with given prefix and suffix.
#
macro(isaac_find_any_library name header library prefix suffix)
if    (NOT ${name}_INCLUDE_DIR)
    CHECK_INCLUDE_FILE(${header} ${name}_INCLUDE_DIR)
endif (NOT ${name}_INCLUDE_DIR)
if    (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
    message (STATUS "Checking library: ${prefix}${library}${suffix}")
    find_library(${name}_LIBRARY 
                 NAMES "${prefix}${library}${suffix}"
                 HINTS ENV LIBRARY_PATH)
endif (${name}_INCLUDE_DIR AND NOT ${name}_LIBRARY)
if(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
    set (HAVE_${name} ${${name}_LIBRARY})
    message (STATUS "Found ${name}  header: ${${name}_INCLUDE_DIR}/${header}")
    message (STATUS "Found ${name} library: ${${name}_LIBRARY}")
endif(${name}_INCLUDE_DIR AND ${name}_LIBRARY)
endmacro(isaac_find_any_library)

#   
# Find static or dynamic library according to globally set preferences
#
macro(isaac_find_library name header library)
    isaac_find_any_library("${name}" "${header}" "${library}" "${iSAAC_LIBRARY_PREFIX}" "${iSAAC_LIBRARY_SUFFIX}")
endmacro(isaac_find_library)

include(CheckIncludeFile)
#   
# Macro to find libraries, with support for static-only search
#
macro(isaac_find_header_or_die variable file)
CHECK_INCLUDE_FILE(${file} ${variable})
if    (${variable})
    message(STATUS "${file} found as ${${variable}}")
else  (${variable})
    message(FATAL_ERROR "Required header ${file} not found.")
endif (${variable})
endmacro(isaac_find_header_or_die)

