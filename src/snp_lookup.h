/**
 * Copyright (c) 2012-2013, Mattias Frånberg
 * All rights reserved.
 *
 * This file is distributed under the Modified BSD License. See the COPYING file
 * for details.
 */

#ifndef __SNP_LOOKUP_H__
#define __SNP_LOOKUP_H__

#ifdef __cplusplus
extern "C" {
#endif

#if HAVE_ENDIAN_H
#include <endian.h>
#elif HAVE_MACHINE_ENDIAN_H
#include <machine/endian.h>
#elif HAVE_SYS_ENDIAN_H
#include <sys/endian.h>
#endif

/**
 * Maps an unpacked snp to its corresponding bits.
 */
unsigned char snp_to_bits[] = { 0, 2, 3, 1 };

/**
 * This files contains a lookup table that maps
 * SNPs packed in a single byte into an array of
 * four bytes.
 */
union snp_lookup_t
{
    /**
     * Accessible as an array.
     */
    unsigned char snp_array[4];

    /**
     * Accessible as a block of bytes.
     */
    int32_t snp_block;
};

#if __BYTE_ORDER == __LITTLE_ENDIAN
#include "snp_lookup_little.h"
#else
#include "snp_lookup_big.h"
#endif /* End test endianess */

#ifdef __cplusplus
}
#endif

/**
 * Integral type used for storing a single SNP.
 */
typedef unsigned char snp_t;

/**
 * Take an unpacked array of SNPs where each SNP is packed in 2 bits,
 * and unpack into a byte array. This function assumes that the bits
 * are packed in the following manner:
 * - Each byte contains 4 SNPs
 * - The SNPs are read from right to left in each byte.
 * - The packed SNPs encoded as follows:
 *   * 00 is homozygous major 
 *   * 01 is missing value
 *   * 10 is hetrozygous
 *   * 11 is homozygous minor
 *
 * - The unpacked SNPs are encoded as follows:
 *   * 0 is homozygous major
 *   * 1 is hetrozygous 
 *   * 2 is homozygous minor
 *   * 3 is missing value
 *
 * @param packed_snps The packed SNPs.
 * @param unpacked_snps The unpacked SNPs.
 * @param num_cols The number of SNPs. 
 */
void
unpack_snps(const snp_t *packed_snps, unsigned char *unpacked_snps, size_t num_cols)
{
    int index;
    int packed_left;

    /* Unpack SNPs in pairs of 4. */
    int32_t *unpacked_snps_p = (int32_t *) unpacked_snps;
    int i;
    int packed_length = num_cols / 4;
    for(i = 0; i < packed_length; i++)
    { 
        *unpacked_snps_p = snp_lookup[ packed_snps[ i ] ].snp_block;
        unpacked_snps_p += 1;
    }

    /* Unpack the trailing SNPs */
    index = packed_length * 4;
    packed_left = num_cols % 4;
    for(i = 0; i < packed_left; i++)
    {
        unpacked_snps[ index + i ] = snp_lookup[ packed_snps[ packed_length ] ].snp_array[ i ];
    }
}

/**
 * Does the reverse of unpack_snps, and packs SNPs into bytes. See unpack_snps for
 * the detailed format.
 *
 * @param unpack_snps The unpacked SNPs.
 * @param packed_snps The packed SNPs.
 * @param num_cols The number of columns.
 */
void
pack_snps(const snp_t *unpacked_snps, unsigned char *packed_snps, size_t num_cols)
{
    int i;
    int packed_index;
    int position_in_byte;

    bzero( packed_snps, (num_cols + 3) / 4 );
    for(i = 0; i < num_cols; i++)
    {
        /* Genotypes are stored backwards. */
        packed_index = i / 4;
        position_in_byte = (i % 4) * 2;
        packed_snps[ packed_index ] |= ( snp_to_bits[ unpacked_snps[ i ] ] << position_in_byte );
    }
}

#endif /* End of __SNP_LOOKUP_H__ */
