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
 ** \file Threads.cpp
 **
 ** Thread management helper utilities.
 **
 ** \author Roman Petrovski
 **/

#include "common/Threads.hpp"

namespace isaac
{
namespace common
{

__thread int runOnNode_ = -1;


} // namespace common
} // namespace isaac
