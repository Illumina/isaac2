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
 ** \file Quality.hh
 **
 ** \brief Various functions and tables to support alignment quality.
 **
 ** \author Come Raczy
 **/

#ifndef iSAAC_ALIGNMENT_QUALITY_HH
#define iSAAC_ALIGNMENT_QUALITY_HH

#include <limits>
#include <vector>

#include <boost/format.hpp>

#include "alignment/Cluster.hh"
#include "common/MathCompatibility.hh"

namespace isaac
{
namespace alignment
{

/**
 ** \brief Utility class providing various services related to base and sequence quality
 **/
class Quality
{
public:
    /**
     ** \brief Return the natural log of the probability of an incorrect base for a
     ** given quality.
     **
     ** This is simply 'log(perror)' where 'perror=10^(-quality/10)'
     **
     ** \param quality PHRED quality score
     **/
    static double getLogError(const unsigned int quality)
    {
        ISAAC_ASSERT_MSG(quality < logErrorLookup.size(),
                         (boost::format("Incorrect quality %u ") % quality).str().c_str());
        return logErrorLookup[quality];
    }

    /**
     ** \brief Return the natural log of the probability of a correct base for a
     ** given quality.
     **
     ** This is simply 'log(1-perror)' where 'perror=10^(-quality/10)'
     **
     ** \param quality PHRED quality score
     **/
    static double getLogCorrect(const unsigned int quality)
    {
        ISAAC_ASSERT_MSG(quality < logMatchLookup.size(),
                         (boost::format("Incorrect quality %u ") % quality).str().c_str());
        return logMatchLookup[quality];
    }

    /**
     ** \brief Return the natural log of the probability of a base that matches reference to be correct.
     **
     ** This is simply 'log(1-perror)' where 'perror=10^(-quality/10)'
     **
     ** \param quality PHRED quality score
     **/
    static double getLogMatch(const unsigned int quality)
    {
        return getLogCorrect(quality);
    }


    /**
     * \brief Same as getLogMismatchSlow but uses a pre-built lookup table
     */
    static double getLogMismatch(const unsigned int quality)
    {
        ISAAC_ASSERT_MSG(quality < logMatchLookup.size(),
                         (boost::format("Incorrect quality %u ") % quality).str().c_str());
        return logMismatchLookup[quality];
    }

    /**
     ** \brief Return the natural log of the probability of a base that mismatches the reference to be wrong.
     ** 
     ** This is 'log(perror/3/pcorrect)' where 'perror=10^(-quality/10)' and
     ** 'pcorrect=1-perror'. The rationale is that if there is an error, each of
     ** the three other bases has 1/3 of the chances of being the correct one.
     ** 
     ** \param quality PHRED quality score
     **/
    static double getLogMismatchSlow(const unsigned quality)
    {
        const double mismatch = pow(10.0, (double)quality / -10.0);
        return log(mismatch / 3.0);
    }

    /**
     ** \brief Return the 'rest of the genome' correction for uniquely aligned reads
     **
     **/
    static double restOfGenomeCorrection(const unsigned genomeLength, const unsigned readLength)
    {
        using std::log;
        return exp(log(2.0) + log((double)genomeLength) - (log(4.0) * (double)readLength));
    }

private:
    /// lookup for log of probability of a sequencing error for a given quality
    static const std::vector<double> logErrorLookup;
    /// lookup for log of probability of a match for a given quality
    static const std::vector<double> logMatchLookup;
    /// lookup for log of probability of a mismatch for a given quality
    static const std::vector<double> logMismatchLookup;
};

void trimLowQualityEnds(Cluster &cluster, const unsigned baseQualityCutoff);


template <typename FpT> bool ISAAC_LP_EQUALS(FpT left, FpT right)
{
    return 0.0000001 >= std::fabs(left - right);
}

template <typename FpT> bool ISAAC_LP_LESS(FpT left, FpT right)
{
    return !ISAAC_LP_EQUALS(left, right) && left < right;
}

} // namespace alignemnt
} // namespace isaac

#endif // #ifndef iSAAC_ALIGNMENT_QUALITY_HH
