<?xml version="1.0"?>
<!--
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
 ** \file GetSupportedSeedLengths.xsl
 **
 ** \author Roman Petrovski
 **/
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
>

<xsl:output method="text"/>

<xsl:template match="/"> 
    <xsl:for-each select="SortedReference/Permutations/Permutation/Masks">
        <xsl:value-of select="@SeedLength"/><xsl:text>
</xsl:text>
    </xsl:for-each>
</xsl:template>

</xsl:stylesheet>
