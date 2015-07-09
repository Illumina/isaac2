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
 ** \file BinLoader.hh
 **
 ** Loads single alignment bin into memory.
 ** 
 ** \author Roman Petrovski
 **/

#ifndef iSAAC_BUILD_BIN_LOADER_HH
#define iSAAC_BUILD_BIN_LOADER_HH

#include <numeric>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/BinMetadata.hh"
#include "build/FragmentIndex.hh"
#include "build/BinData.hh"

namespace isaac
{
namespace build
{

class BinLoader
{
public:
    BinLoader()
    {
    }

    void loadData(BinData &data);

private:
    void loadUnalignedData(BinData &binData);
    void loadAlignedData(BinData &binData);
    const io::FragmentAccessor &loadFragment(BinData &binData, std::istream &isData, unsigned long &offset);
    void storeFragmentIndex(const io::FragmentAccessor& mateFragment,
                            unsigned long mateOffset, unsigned long offset,
                            BinData& binData);
};


} // namespace build
} // namespace isaac

#endif // #ifndef iSAAC_BUILD_BIN_LOADER_HH
