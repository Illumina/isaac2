/**
 ** Isaac Genome Alignment Software
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
 ** \file KUniqueness.hh
 **
 ** \brief Basic declarations for computing k-uniqueness
 **
 ** \author Come Raczy
 **/

#ifndef ISAAC_REFERENCE_K_UNIQUENESS_HH
#define ISAAC_REFERENCE_K_UNIQUENESS_HH

namespace isaac
{
namespace reference
{

typedef unsigned short AnnotationValue;
typedef std::vector<AnnotationValue> Annotation;
typedef std::vector<AnnotationValue> ContigAnnotation;
typedef std::vector<ContigAnnotation> ContigAnnotations;
typedef std::vector<ContigAnnotations> ContigAnnotationsList;

typedef AnnotationValue DistanceToBeNeighborless;
static const DistanceToBeNeighborless K_UNIQUE_TOO_FAR = DistanceToBeNeighborless(0) - 1;

} // namespace reference
} //namespace isaac

#endif // #ifndef ISAAC_REFERENCE_K_UNIQUENESS_HH
