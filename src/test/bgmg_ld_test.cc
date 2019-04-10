#include "gtest/gtest.h"
#include "omp.h"

#include <cstdlib>
#include <set>
#include <iostream>
#include <random>
#include <algorithm>
#include <string>

#include "bgmg_parse.h"
#include "plink_ld.h"
#include "snp_lookup.h"

const std::string DataFolder = "/home/oleksanf/github/mixer/src/testdata";

// interesting detail: in LD r2 calculation plink calculates the mean across
// genotypes defined in both SNPs --- therefore we pass the "second" argument below. 
double mean(const std::string& unpacked, const std::string& second) {
  // 0 is homozygous major, 1 is hetrozygous,  2 is homozygous minor, 3 is missing value
  double sum = 0.0, count = 0.0;
  for (int i = 0; i < unpacked.size(); i++) {
    if (unpacked[i] == 3 || second[i] == 3) continue;
    sum += unpacked[i];
    count += 1;
  }
  return sum / count;
}

double corr(const std::string& unpacked1, const std::string& unpacked2) {
  // 0 is homozygous major, 1 is hetrozygous,  2 is homozygous minor, 3 is missing value
  double mean1 = mean(unpacked1, unpacked2);
  double mean2 = mean(unpacked2, unpacked1);
  double cov = 0.0, var1 = 0.0, var2 = 0.0;
  for (int i = 0; i < unpacked1.size(); i++) {
    if (unpacked1[i] == 3 || unpacked2[i] == 3) continue;
    double diff1 = unpacked1[i] - mean1;
    double diff2 = unpacked2[i] - mean2;
    cov += diff1 * diff2;
    var1 += diff1 * diff1;
    var2 += diff2 * diff2;
  }
  return cov / sqrt(var1 * var2);
}

void generate_genotypes(int num_subj, int num_snps, double missing_rate, std::vector<std::string>*unpacked_snps, std::string* buffer) {
  int unfiltered_sample_ct4 = (num_subj + 3) / 4;

  for (int i = 0; i < num_snps; i++) {
    std::string unpacked_snp(num_subj, (char)0);
    for (int j = 0; j < num_subj; j++) {
      if (((float)rand() / (float)RAND_MAX) < missing_rate)
        unpacked_snp[j] = 3;
      else
        unpacked_snp[j] = (char)(rand() % 3);
    }
    unpacked_snps->push_back(unpacked_snp);
 }

  buffer->assign(3, char(0));
  std::string packed_snp(unfiltered_sample_ct4, (char)0);
  std::string unpacked_snp(num_subj, (char)0);
  for (int i = 0; i < num_snps; i++) {
    pack_snps((const snp_t*)unpacked_snps->at(i).c_str(), (unsigned char*)&packed_snp[0], num_subj);
    buffer->append(packed_snp);

    // tests that pack&unpack gives the original genotype vector
    unpack_snps((const snp_t *)&packed_snp[0], (unsigned char*)&unpacked_snp[0], num_subj);
    for (int j = 0; j < num_subj; j++)
      ASSERT_EQ(unpacked_snps->at(i)[j], unpacked_snp[j]);
  }

  // std::string buffer2(3 + unfiltered_sample_ct4 * num_snps, 0);
  // for (int i = 0; i < buffer2.size(); i++) buffer2[i] = rand() % 256;
}

void test_ld(int num_subj, int num_snps, double missing_rate) {
  std::vector<std::string> unpacked_snps;
  std::string buffer;
  generate_genotypes(num_subj, num_snps, missing_rate, &unpacked_snps, &buffer);

  FILE* bedfile = fmemopen(&buffer[0], buffer.size(), "rb");
  PlinkLdBedFileChunk plink_ld(num_subj, 0, num_snps, bedfile);
  for (int i = 0; i < num_snps; i++) {
    for (int j = i; j < num_snps; j++) {
      double rIJ = PlinkLdBedFileChunk::calculate_ld_corr(plink_ld, plink_ld, i, j);
      double rJI = PlinkLdBedFileChunk::calculate_ld_corr(plink_ld, plink_ld, j, i);
      double rIJ_plain = corr(unpacked_snps[i], unpacked_snps[j]);
      if (std::isfinite(rIJ) || std::isfinite(rJI) || std::isfinite(rIJ_plain)) {
        ASSERT_EQ(rIJ, rJI);
        ASSERT_NEAR(rIJ, rIJ_plain, 1e-12);
      }
    }
  }
  fclose(bedfile);
}

// --gtest_filter=TestLd.*
TEST(TestLd, WithoutMissingGenotypes) {
  for (int num_subj = 1; num_subj < 100; num_subj++) test_ld(num_subj, 10, 0.0);
  for (int num_subj = 10001; num_subj < 10100; num_subj++) test_ld(num_subj, 10, 0.0);
}

TEST(TestLd, RandomMissingGenotypes) {
  for (int num_subj = 1; num_subj < 100; num_subj++) test_ld(num_subj, 10, 0.5);
  for (int num_subj = 10001; num_subj < 10100; num_subj++) test_ld(num_subj, 10, 0.5);
}

TEST(TestLd, AllMissingGenotypes) {
  for (int num_subj = 1; num_subj < 100; num_subj++) test_ld(num_subj, 10, 1.001);
  for (int num_subj = 10001; num_subj < 10100; num_subj++) test_ld(num_subj, 10, 1.001);
}

TEST(TestLd, DifferentChunks) {
  std::vector<std::string> unpacked_snps;
  std::string buffer;
  int num_subj = 789, num_snps=20;
  generate_genotypes(num_subj, num_snps, 0.1, &unpacked_snps, &buffer);
  FILE* bedfile = fmemopen(&buffer[0], buffer.size(), "rb");

  PlinkLdBedFileChunk big_chunk(num_subj, 0, num_snps, bedfile);
  PlinkLdBedFileChunk left_chunk(num_subj, 0, num_snps/2, bedfile);
  PlinkLdBedFileChunk right_chunk(num_subj, num_snps/2, num_snps/2, bedfile);

  for (int i = 0; i < num_snps/2; i++) {
    for (int j = 0; j < num_snps/2; j++) {
      double rIJ_big = PlinkLdBedFileChunk::calculate_ld_corr(big_chunk, big_chunk, i, num_snps/2 + j);
      double rIJ_left = PlinkLdBedFileChunk::calculate_ld_corr(left_chunk, right_chunk, i, j);
      ASSERT_TRUE(std::isfinite(rIJ_big) && std::isfinite(rIJ_left));
      ASSERT_EQ(rIJ_big, rIJ_left);
    }
  }

  fclose(bedfile);
}
