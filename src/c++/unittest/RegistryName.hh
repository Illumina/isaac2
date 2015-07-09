/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 ** \file RegistryName.hh
 **
 ** Management of the registry names for the cppunit tests.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_UNIT_TEST_REGISTRY_NAME
#define iSAAC_UNIT_TEST_REGISTRY_NAME

#include <stdexcept>
#include <string>
#include <vector>

const std::vector<std::string> &getRegistryNameList();
std::string registryName(const std::string &name) throw (std::invalid_argument);

#endif // #ifndef iSAAC_UNIT_TEST_REGISTRY_NAME
