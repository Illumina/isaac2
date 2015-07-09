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
 ** \file RegistryName.cpp
 **
 ** Management of the registry names for the cppunit tests.
 **
 ** \author Come Raczy
 **/

#include "RegistryName.hh"

#include <boost/filesystem.hpp>
#include <vector>
#include <string>
#include <fstream>

boost::filesystem::path getFilePath() {return "RegistryNames.txt";}

std::vector<std::string> initializeNameList()
{
    std::vector<std::string> nameList;
    if(boost::filesystem::exists(getFilePath()))
    {
        std::ifstream is(getFilePath().string().c_str());
        std::string name;
        while (getline(is, name))
        {
            if (!name.empty() && nameList.end() == std::find(nameList.begin(), nameList.end(), name))
            {
                nameList.push_back(name);
            }
        }
    }
    return nameList;
}

const std::vector<std::string> &getRegistryNameList()
{
    static const std::vector<std::string> nameList = initializeNameList();
    return nameList;
}

std::string registryName(const std::string &name) throw (std::invalid_argument)
{
    const std::vector<std::string> nameList = getRegistryNameList();
    const std::vector<std::string>::const_iterator found = std::find(nameList.begin(), nameList.end(), name);
    if (found != nameList.end())
        return name;

    const std::vector<std::string>::const_iterator commentedFound = std::find(nameList.begin(), nameList.end(), "#" + name);
    if (commentedFound != nameList.end())
        return name;
    else
        throw std::invalid_argument(std::string("Not a registryName: ") + name +
                                    std::string(" [check that ") + getFilePath().string() +
                                    std::string(" constains '") + name + ". Use #"+ name + " to temporary disable a test" +
                                    std::string("']"));
}
