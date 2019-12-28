/**
 * Copyright (c) 2012-2013, Mattias Frånberg
 * All rights reserved.
 *
 * This file is distributed under the Modified BSD License. See the COPYING file
 * for details.
 */

#ifndef __SNP_LOOKUP_H__
#define __SNP_LOOKUP_H__

#include "stdlib.h"  // size_t

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

#if __BYTE_ORDER != __LITTLE_ENDIAN
#error Big endian architectures are not supported
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
unpack_snps(const snp_t *packed_snps, unsigned char *unpacked_snps, size_t num_cols);

/**
 * Does the reverse of unpack_snps, and packs SNPs into bytes. See unpack_snps for
 * the detailed format.
 *
 * @param unpack_snps The unpacked SNPs.
 * @param packed_snps The packed SNPs.
 * @param num_cols The number of columns.
 */
void
pack_snps(const snp_t *unpacked_snps, unsigned char *packed_snps, size_t num_cols);

#endif /* End of __SNP_LOOKUP_H__ */
