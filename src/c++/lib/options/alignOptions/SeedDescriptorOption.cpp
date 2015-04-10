/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Public License 1
 **
 ** You should have received a copy of the Illumina Public License 1
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file SeedDescriptorOption.cpp
 **
 ** seeds option parsing
 **
 ** \author Roman Petrovski
 **/
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "alignment/SeedMetadata.hh"
#include "common/Debug.hh"
#include "common/Exceptions.hh"
#include "oligo/Kmer.hh"

#include "SeedDescriptorOption.hh"


namespace isaac
{
namespace options
{
namespace alignOptions
{

unsigned parseManualSeedDescriptor(
    const std::string &descriptor,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    unsigned ret = 0;
    using boost::algorithm::split;
    using boost::algorithm::is_any_of;

    std::vector<std::string> offsetList;
    split(offsetList, descriptor,  is_any_of(":"));
    if (offsetList.empty())
    {
        const boost::format message = boost::format(
            "\n   *** The list of seed offsets for %s is empty. At least one seed is needed for each read ***\n") %
            readMetadata;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    BOOST_FOREACH(const std::string &offsetString, offsetList)
    {
        try
        {
            const unsigned offset = boost::lexical_cast<unsigned>(offsetString);
            alignment::SeedMetadata seedMetadata(offset, seedLength, readMetadata.getIndex(), seedMetadataList.size());
            if (offset + seedLength > readMetadata.getLength())
            {
                ISAAC_THREAD_CERR << "WARNING: ignored " << seedMetadata <<
                    " as it stretches beyond the read " << readMetadata.getNumber() <<
                    " which is " << readMetadata.getLength() << " bases long" << std::endl;
            }
            else
            {
                seedMetadataList.push_back(seedMetadata);
                ++ret;
                ISAAC_THREAD_CERR << "constructed " << seedMetadataList.back() << std::endl;
            }
        }
        catch(boost::bad_lexical_cast &)
        {
            const boost::format message = boost::format("\n   *** Invalid seed offset '%s' found in '%s' ***\n") %
                offsetString % descriptor;
            BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
        }
    }
    return ret;
}


/**
 * \return number of seeds generated
 */
unsigned parseAutoSeedDescriptor(
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    unsigned generated = 0;
    unsigned offset = 0;
    unsigned endOffset = readMetadata.getLength();

    //put two seeds at the extremities so that we get best chance of missing the homopolymers
    if (readMetadata.getLength() > seedLength)
    {
        {
            alignment::SeedMetadata seedMetadata(0, seedLength, readMetadata.getIndex(), seedMetadataList.size());
            offset = seedLength;
            seedMetadataList.push_back(seedMetadata);
        }
        ISAAC_THREAD_CERR << "constructed extremity seed " << seedMetadataList.back() << std::endl;

        {
            endOffset = readMetadata.getLength()  - seedLength;
            alignment::SeedMetadata seedMetadata(endOffset, seedLength, readMetadata.getIndex(), seedMetadataList.size());
            seedMetadataList.push_back(seedMetadata);
        }
        ISAAC_THREAD_CERR << "constructed extremity seed " << seedMetadataList.back() << std::endl;

        generated = 2;
    }

    //put as many seeds as we can in what's left. Don't put spaces between seeds mainly because it is easier to debug this way
    while(offset + seedLength <= endOffset)
    {
        alignment::SeedMetadata seedMetadata(offset, seedLength, readMetadata.getIndex(), seedMetadataList.size());
        seedMetadataList.push_back(seedMetadata);
        ++generated;
        offset += seedLength;
        ISAAC_THREAD_CERR << "constructed " << seedMetadataList.back() << std::endl;
    }

    // if less than 4 seeds generated, add overlapping ones
    offset = seedLength / 2;
    if (endOffset > seedLength / 2)
    {
        endOffset -= seedLength / 2;
        while (generated < 4 && offset + seedLength <= endOffset)
        {
            alignment::SeedMetadata seedMetadata(offset, seedLength, readMetadata.getIndex(), seedMetadataList.size());
            seedMetadataList.push_back(seedMetadata);
            ++generated;
            offset += seedLength;
            ISAAC_THREAD_CERR << "constructed overlapping " << seedMetadataList.back() << std::endl;
        }
    }
    return generated;
}


/**
 * \return Number of generated seeds
 */
unsigned parseStepSeedDescriptor(
    const unsigned step,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    ISAAC_ASSERT_MSG(readMetadata.getLength() >= seedLength, "Read is too short for seed lenght " << seedLength << " " << readMetadata);

    unsigned i = 0;
    for (; i < readMetadata.getLength() - seedLength; i += step)
    {
        alignment::SeedMetadata seedMetadata(i, seedLength, readMetadata.getIndex(), seedMetadataList.size());
        seedMetadataList.push_back(seedMetadata);
    }
    return i;
}

/**
 * \return Maximum number of first pass seeds possible
 */
unsigned parseReadSeedDescriptor(
    const std::string &descriptor,
    const flowcell::ReadMetadata& readMetadata,
    const unsigned seedLength,
    alignment::SeedMetadataList &seedMetadataList)
{
    if (descriptor.empty())
    {
        const boost::format message = boost::format(
            "\n   *** The seed descriptor for %s is empty. At least one seed is needed ***\n") %
            readMetadata;
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    if ("step=" == descriptor.substr(0, 5))
    {
        return parseStepSeedDescriptor(boost::lexical_cast<unsigned>(descriptor.substr(5)), readMetadata, seedLength, seedMetadataList);
    }
    else if ("all" == descriptor)
    {
        return parseStepSeedDescriptor(1, readMetadata, seedLength, seedMetadataList);
    }
    else if ("auto" == descriptor)
    {
        return parseAutoSeedDescriptor(readMetadata, seedLength, seedMetadataList);
    }
    else
    {
        return parseManualSeedDescriptor(descriptor, readMetadata, seedLength, seedMetadataList);
    }
}

/**
 * Parses seed descriptor into seed objects. Reduces firstPassSeeds if the read length does not permit.
 */
alignment::SeedMetadataList parseSeedDescriptor(
    const std::vector<flowcell::ReadMetadata> &readMetadataList,
    const std::string &seedDescriptor,
    const unsigned seedLength,
    unsigned &firstPassSeeds)
{
    if (seedDescriptor.empty())
    {
        const boost::format message = boost::format("\n   *** The seed descriptor is empty. At least one seed is needed ***\n");
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    alignment::SeedMetadataList seedMetadataList;
    std::vector<std::string> seedDescriptorList; // split by read
    using boost::algorithm::split;
    using boost::algorithm::is_any_of;
    split(seedDescriptorList, seedDescriptor,  is_any_of(","));
    if (readMetadataList.size() < seedDescriptorList.size())
    {
        const boost::format message = boost::format("\n   *** Too many lists-of-seeds in seed-descriptor '%s': found %d: %d reads only ***\n") %
            seedDescriptor % seedDescriptorList.size() % readMetadataList.size();
        BOOST_THROW_EXCEPTION(common::InvalidOptionException(message.str()));
    }
    // extend the last list-of-seeds to all subsequent reads if needed
    seedDescriptorList.resize(readMetadataList.size(), seedDescriptorList.back());
    // create all the seeds
    std::vector<flowcell::ReadMetadata>::const_iterator readMetadata = readMetadataList.begin();
    BOOST_FOREACH(const std::string &descriptor, seedDescriptorList)
    {
        const unsigned parsedSeeds = parseReadSeedDescriptor(descriptor, *readMetadata, seedLength, seedMetadataList);
        if (parsedSeeds < firstPassSeeds)
        {
            firstPassSeeds = parsedSeeds;
            ISAAC_THREAD_CERR << "WARNING: reducing first pass seeds to " << firstPassSeeds << std::endl;
        }
        ++readMetadata;
    }

    return seedMetadataList;
}

} // namespace alignOptions
} // namespace option
} // namespace isaac
