#include "gtest/gtest.h"
#include "omp.h"

#include <iostream>
#include <random>
#include <algorithm>

#include "bgmg_calculator.h"

/*
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
*/

namespace {

TEST(BgmgMainTest, ShouldSucceed) {
}

/*
TEST(BgmgGzipTest, TestGzip) {
  std::ifstream file("hello.z", std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::zlib_decompressor());
  in.push(file);
  boost::iostreams::copy(in, std::cout);
}
*/

TEST(BgmgTest, LoadData) {
  // Input data:
  // - LD structure as a set of pairwise LD r2 values
  //   i1, i2, r2 <- currently this came from plink format (e.i. lower triangular, not limited to tag SNPs)
  //   Split per chromosome
  // - LD structure in a most convenient ways
  // - Plink LD structure without differentiating tag and SNP variants
  // - Several vectors, potentially with different length (#tag, #snp)
  //   zvec, nvec, hvec, weights
  // 
}

class TestMother {
public:
  TestMother(int num_snp, int num_tag, int n) : num_snp_(num_snp), num_tag_(num_tag), rd_(), g_(rd_()) {
    std::uniform_real_distribution<> hvec_dist(0.0f, 0.5f); // 2*p(1-p).

    for (int i = 0; i < num_snp_; i++) tag_to_snp_.push_back(i);
    std::shuffle(tag_to_snp_.begin(), tag_to_snp_.end(), g_);
    tag_to_snp_.resize(num_tag_);

    regenerate_zvec();
    for (int i = 0; i < num_tag_; i++) n_vec_.push_back(n);
    for (int i = 0; i < num_tag_; i++) weights_.push_back(1.0f);
    for (int i = 0; i < num_snp_; i++) h_vec_.push_back(hvec_dist(g_));
  }

  std::vector<int>* tag_to_snp() { return &tag_to_snp_; }
  std::vector<float>* zvec() { return &z_vec_; }
  std::vector<float>* nvec() { return &n_vec_; }
  std::vector<float>* hvec() { return &h_vec_; }
  std::vector<float>* weights() { return &weights_; }

  void regenerate_zvec() {
    z_vec_.clear();
    std::normal_distribution<float> norm_dist(0.0, 1.5);
    for (int i = 0; i < num_tag_; i++) z_vec_.push_back(norm_dist(g_));
  }

  void make_r2(int num_r2, std::vector<int>* snp_index, std::vector<int>* tag_index, std::vector<float>* r2) {
    std::uniform_int_distribution<> dis(0, num_snp_ - 1);
    std::uniform_real_distribution<> dis2(0.0f, 1.0f);
    while (tag_index->size() < num_r2) {
      int tag = dis(g_);
      int snp = dis(g_);
      if (tag <= snp) continue;
      tag_index->push_back(tag);
      snp_index->push_back(snp);
      r2->push_back(dis2(g_));
    }
  }

  std::mt19937& random_engine() { return g_; }
private:
  int num_snp_;
  int num_tag_;
  std::vector<int> tag_to_snp_;
  std::vector<float> z_vec_;
  std::vector<float> h_vec_;
  std::vector<float> n_vec_;
  std::vector<float> weights_;
  std::random_device rd_;
  std::mt19937 g_;
};

TEST(UgmgTest, CalcLikelihood) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);

  int trait = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));
  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);
  
  calc.set_ld_r2_coo(r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
                  
  // TBD: validate CSR structure (set_ld_r2_csr is quite tricky)

  calc.set_hvec(num_snp, &tm.hvec()->at(0));

  if (true) {
    //calc.find_snp_order();

    // Now, let's calculate log likelihood.
    calc.find_tag_r2sum(0, 2.1);
    calc.find_tag_r2sum(0, 1.1);
    calc.find_tag_r2sum(0, 3.1);
    calc.find_tag_r2sum(0, 3);
    calc.find_tag_r2sum(0, 4);
    calc.find_tag_r2sum(0, 3.9);
  }

  double cost = calc.calc_univariate_cost(0.2, 1.2, 0.1);
  double cost_nocache = calc.calc_univariate_cost_nocache(0.2, 1.2, 0.1);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_FLOAT_EQ(cost, cost_nocache);

  std::vector<float> zvec_grid, zvec_pdf, zvec_pdf_nocache;
  for (float z = 0; z < 15; z += 0.1) {
    zvec_grid.push_back(z);
    zvec_pdf.push_back(0.0f);
    zvec_pdf_nocache.push_back(0.0f);
  }

  calc.calc_univariate_pdf(0.2, 1.2, 0.1, zvec_grid.size(), &zvec_grid[0], &zvec_pdf[0]);
  calc.set_option("diag", 0.0);

  calc.set_option("cache_tag_r2sum", 0);
  calc.calc_univariate_pdf(0.2, 1.2, 0.1, zvec_grid.size(), &zvec_grid[0], &zvec_pdf_nocache[0]);

  for (int i = 0; i < zvec_pdf_nocache.size(); i++)
    ASSERT_FLOAT_EQ(zvec_pdf[i], zvec_pdf_nocache[i]);

  calc.set_option("fast_cost", 1);
  cost = calc.calc_univariate_cost(0.2, 1.2, 0.1);
  ASSERT_TRUE(std::isfinite(cost));
}

TEST(BgmgTest, CalcLikelihood) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 3);
  calc.set_option("cache_tag_r2sum", 1);

  int trait = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  trait = 2; tm.regenerate_zvec();
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);

  calc.set_ld_r2_coo(r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

  calc.set_weights_randprune(20, 0.25);
  calc.set_hvec(num_snp, &tm.hvec()->at(0));

  float pi_vec[] = { 0.1, 0.2, 0.15 };
  float sig2_beta[] = { 0.5, 0.3 };
  float rho_beta = 0.8;
  float sig2_zero[] = { 1.1, 1.2 };
  float rho_zero = 0.1;

  double cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  double cost_nocache = calc.calc_bivariate_cost_nocache(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_FLOAT_EQ(cost, cost_nocache);

  calc.set_option("diag", 0.0);

  calc.set_option("fast_cost", 1);
  cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  ASSERT_TRUE(std::isfinite(cost));

  std::vector<float> zvec1_grid, zvec2_grid, zvec_pdf, zvec_pdf_nocache;
  for (float z1 = -10; z1 < 10; z1 += 0.2) {
    for (float z2 = -10; z2 < 10; z2 += 0.2) {
      zvec1_grid.push_back(z1);
      zvec2_grid.push_back(z2);
      zvec_pdf.push_back(0.0f);
      zvec_pdf_nocache.push_back(0.0f);
    }
  }

  calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, zvec_pdf.size(), &zvec1_grid[0], &zvec2_grid[0], &zvec_pdf[0]);

  calc.set_option("cache_tag_r2sum", 0);
  calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, zvec_pdf.size(), &zvec1_grid[0], &zvec2_grid[0], &zvec_pdf_nocache[0]);

  for (int i = 0; i < zvec_pdf_nocache.size(); i++)
    ASSERT_FLOAT_EQ(zvec_pdf[i], zvec_pdf_nocache[i]);
}

// bgmg-test.exe --gtest_filter=Test.RandomSeedAndThreading
TEST(Test, RandomSeedAndThreading) {
  int num_snp = 100;
  int num_tag = 50;
  int kmax = 200; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int num_r2 = 20;
  TestMother tm(num_snp, num_tag, N);
  TestMother tm2(num_snp, num_tag, N);
  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(num_r2, &snp_index, &tag_index, &r2);

  int64_t seed_list[3] = { 123123123, 456456456, 123123123 };

  double ugmg_costs[3];
  double bgmg_costs[3];
  for (int num_threads = 1; num_threads <= 16; num_threads *= 2) {
    for (int i = 0; i < 3; i++) {
      BgmgCalculator calc;
      calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
      calc.set_option("max_causals", num_snp);
      calc.set_option("kmax", kmax);
      calc.set_option("num_components", 3);
      calc.set_option("cache_tag_r2sum", 1);
      calc.set_option("threads", num_threads);
      calc.set_seed(seed_list[i]);

      int num_threads_real;
#pragma omp parallel
{
      num_threads_real = omp_get_num_threads();
}
      ASSERT_EQ(num_threads_real, num_threads);

      int trait = 1;
      calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
      calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));
      trait = 2;
      calc.set_zvec(trait, num_tag, &tm2.zvec()->at(0));
      calc.set_nvec(trait, num_tag, &tm2.nvec()->at(0));

      calc.set_ld_r2_coo(r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
      calc.set_ld_r2_csr();  // finalize csr structure

      calc.set_weights_randprune(20, 0.25);
      calc.set_hvec(num_snp, &tm.hvec()->at(0));

      float pi_vec[] = { 0.1, 0.2, 0.15 };
      float sig2_beta[] = { 0.5, 0.3 };
      float rho_beta = 0.8;
      float sig2_zero[] = { 1.1, 1.2 };
      float rho_zero = 0.1;
      double ugmg_cost = calc.calc_univariate_cost(pi_vec[0], sig2_zero[0], sig2_beta[0]);
      if (num_threads == 1) ugmg_costs[i] = ugmg_cost;
      else ASSERT_FLOAT_EQ(ugmg_costs[i], ugmg_cost);
      ASSERT_FLOAT_EQ(ugmg_cost, calc.calc_univariate_cost(pi_vec[0], sig2_zero[0], sig2_beta[0]));  // check calc twice => the same cost
      ASSERT_FLOAT_EQ(ugmg_cost, calc.calc_univariate_cost_nocache(pi_vec[0], sig2_zero[0], sig2_beta[0]));

      double bgmg_cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
      if (num_threads == 1) bgmg_costs[i] = bgmg_cost;
      else ASSERT_FLOAT_EQ(bgmg_costs[i], bgmg_cost);
      ASSERT_FLOAT_EQ(bgmg_costs[i], calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero));  // check calc twice => the same cost
      ASSERT_FLOAT_EQ(bgmg_costs[i], calc.calc_bivariate_cost_nocache(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero));
    }
    ASSERT_FLOAT_EQ(ugmg_costs[0], ugmg_costs[2]);
    ASSERT_TRUE(abs(ugmg_costs[0] - ugmg_costs[1]) > 1e-4);
    ASSERT_FLOAT_EQ(bgmg_costs[0], bgmg_costs[2]);
    ASSERT_TRUE(abs(bgmg_costs[0] - bgmg_costs[1]) > 1e-4);
  }
}

// bgmg-test.exe --gtest_filter=Test.tag_r2_caching
TEST(Test, tag_r2_caching) {
  int num_snp = 100;
  int num_tag = 50;
  int kmax = 200; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int num_r2 = 20;
  TestMother tm(num_snp, num_tag, N);
  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(num_r2, &snp_index, &tag_index, &r2);
  
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_seed(123123123);
  int trait = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  int sequence_length = 10;
  std::vector<float> num_causal_sequence;
  std::uniform_real_distribution<float> rng(51.0, 99.5);
  for (int i = 0; i < sequence_length; i++) num_causal_sequence.push_back(rng(tm.random_engine()));
  
  calc.set_ld_r2_coo(r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

  calc.set_weights_randprune(20, 0.25);
  calc.set_hvec(num_snp, &tm.hvec()->at(0));

  std::vector<double> costs(sequence_length, 0.0);
  for (int j = 0; j < 100; j++) {  // repeat the sequence 100 times and validate that we got the same cost.
    for (int i = 0; i < sequence_length; i++) {
      float pi = static_cast<float>(num_causal_sequence[i]) / static_cast<float>(num_snp);
      double cost = calc.calc_univariate_cost(pi, 0.2, 0.15);
      if (j == 0) costs[i] = cost;
      ASSERT_FLOAT_EQ(cost, costs[i]);
    }
  }
}

}  // namespace
