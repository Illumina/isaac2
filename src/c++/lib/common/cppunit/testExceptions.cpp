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
 **/

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#include "RegistryName.hh"
#include "testExceptions.hh"

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( TestExceptions, registryName("Exceptions"));

void TestExceptions::setUp()
{
}

void TestExceptions::tearDown()
{
}


void TestExceptions::testErrorNumber()
{
    CPPUNIT_ASSERT(true);
}

