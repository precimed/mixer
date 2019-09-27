#include "gtest/gtest.h"
#include "omp.h"

#include <set>
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
  TestMother(int num_snp, int num_tag, int n) : num_snp_(num_snp), num_tag_(num_tag), rd_(), g_(12341234) {
    std::uniform_real_distribution<> mafvec_dist(0.0f, 0.5f);

    for (int i = 0; i < num_snp_; i++) tag_to_snp_.push_back(i);
    std::shuffle(tag_to_snp_.begin(), tag_to_snp_.end(), g_);
    tag_to_snp_.resize(num_tag_);
    std::sort(tag_to_snp_.begin(), tag_to_snp_.end()); 

    regenerate_zvec();
    for (int i = 0; i < num_tag_; i++) n_vec_.push_back(n);
    for (int i = 0; i < num_tag_; i++) weights_.push_back(1.0f);
    for (int i = 0; i < num_snp_; i++) mafvec_.push_back(mafvec_dist(g_));
    for (int i = 0; i < num_snp_; i++) chrnumvec_.push_back(1);
  }

  std::vector<int>* tag_to_snp() { return &tag_to_snp_; }
  std::vector<float>* zvec() { return &z_vec_; }
  std::vector<float>* nvec() { return &n_vec_; }
  std::vector<float>* mafvec() { return &mafvec_; }
  std::vector<int>* chrnumvec() { return &chrnumvec_; }
  std::vector<float>* weights() { return &weights_; }

  void regenerate_zvec() {
    z_vec_.clear();
    std::normal_distribution<float> norm_dist(0.0, 1.5);
    for (int i = 0; i < num_tag_; i++) z_vec_.push_back(norm_dist(g_));
  }

  void make_r2(int num_r2, std::vector<int>* snp_index, std::vector<int>* tag_index, std::vector<float>* r2) {
    std::uniform_real_distribution<> dis2(0.0f, 1.0f);
    std::vector<float> rand_list;
    for (int i = 0; i < 1000; i++) rand_list.push_back(((i%3==0) ? -1.0f : 1.0f) * dis2(g_));
    
    std::set<std::tuple<int, int>> pairs;

    std::uniform_int_distribution<> dis(0, num_snp_ - 1);
    int num_try = 0;
    while ((tag_index->size() < num_r2) && (num_try < 10*num_r2)) {
      num_try++;
      int tag = dis(g_);
      int snp = dis(g_);
      if (tag <= snp) continue;
      if (pairs.find(std::make_tuple(tag, snp)) != pairs.end()) continue;
      tag_index->push_back(tag);
      snp_index->push_back(snp);
      pairs.insert(std::make_tuple(tag, snp));
      r2->push_back(rand_list[r2->size() % rand_list.size()]);
    }
  }

  std::mt19937& random_engine() { return g_; }
private:
  int num_snp_;
  int num_tag_;
  std::vector<int> tag_to_snp_;
  std::vector<float> z_vec_;
  std::vector<float> mafvec_;
  std::vector<int> chrnumvec_;
  std::vector<float> n_vec_;
  std::vector<float> weights_;
  std::random_device rd_;
  std::mt19937 g_;
};

// --gtest_filter=LdTest.ValidateMultipleChromosomes
TEST(LdTest, ValidateMultipleChromosomes) {
  int num_snp = 60;
  int num_tag = 40;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(60, 40, N);
  
  std::vector<int> chrnumvec; 
  for (int i = 0; i < 40; i++) chrnumvec.push_back(1);
  for (int i = 0; i < 20; i++) chrnumvec.push_back(2);
  
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("r2min", 0.15f);

  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &chrnumvec[0]);

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(200, &snp_index, &tag_index, &r2);
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);

  calc.set_ld_r2_csr();  // finalize csr structure

  // retrieve cached and non-cached LDr2 sum, and compare the result
  calc.set_option("cache_tag_r2sum", 0);
  std::vector<float> tag_r2_sum(num_tag*kmax, 0.0f);
  calc.retrieve_tag_r2_sum(0, 1, num_tag*kmax, &tag_r2_sum[0]);
  calc.set_option("cache_tag_r2sum", 1);
  std::vector<float> tag_r2_sum_cached(num_tag*kmax, 0.0f);
  calc.retrieve_tag_r2_sum(0, 1, num_tag*kmax, &tag_r2_sum_cached[0]);
  for (int i = 0; i < num_tag*kmax; i++) ASSERT_FLOAT_EQ(tag_r2_sum[i], tag_r2_sum_cached[i]);
}

void UgmgTest_CalcLikelihood(float r2min, int trait_index) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 5;
  int kmax = 20; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);

  calc.set_zvec(trait_index, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait_index, num_tag, &tm.nvec()->at(0));
  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

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

  double cost = calc.calc_univariate_cost(trait_index, 0.2, 1.2, 0.1);
  double cost_nocache = calc.calc_univariate_cost_nocache(trait_index, 0.2, 1.2, 0.1);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_FLOAT_EQ(cost, cost_nocache);

  std::vector<float> zvec_grid, zvec_pdf, zvec_pdf_nocache, zvec_pdf_unified;
  for (float z = 0; z < 15; z += 0.1) {
    zvec_grid.push_back(z);
    zvec_pdf.push_back(0.0f);
    zvec_pdf_nocache.push_back(0.0f);
    zvec_pdf_unified.push_back(0.0f);
  }

  const float sig2_beta = 0.1f;
  const float pi_val = 0.2f;
  const float sig2_zeroL = pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;
  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);

  calc.calc_univariate_pdf(trait_index, pi_val, sig2_zeroA, sig2_beta, zvec_grid.size(), &zvec_grid[0], &zvec_pdf[0]);
  calc.set_option("diag", 0.0);

  calc.set_option("cache_tag_r2sum", 0);
  calc.calc_univariate_pdf(trait_index, pi_val, sig2_zeroA, sig2_beta, zvec_grid.size(), &zvec_grid[0], &zvec_pdf_nocache[0]);

  for (int i = 0; i < zvec_pdf_nocache.size(); i++) {
    ASSERT_NEAR(zvec_pdf[i], zvec_pdf_nocache[i], 2e-7);  // 4.93722e-05 vs 4.9372218e-05 due to approximation of float as uint16_t
  }

  calc.calc_unified_univariate_pdf(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL, zvec_grid.size(), &zvec_grid[0], &zvec_pdf_unified[0]);
  for (int i = 0; i < zvec_pdf_unified.size(); i++) {
    ASSERT_TRUE(std::isfinite(zvec_pdf_unified[i]));
    //  std::cout << zvec_pdf_nocache[i] << ", " << zvec_pdf_unified[i] << "\n";
  }

  //int64_t BgmgCalculator::calc_univariate_power(int trait_index, float pi_vec, float sig2_zero, float sig2_beta, float zthresh, int length, float* nvec, float* svec) {
  std::vector<float> nvec;
  for (int n = 10; n < 1000; n += 10) nvec.push_back(n);
  std::vector<float> svec(nvec.size(), 0.0f);
  float zthresh = 5.45f;

  calc.calc_univariate_power(trait_index, 0.2, 1.2, 0.1, zthresh, nvec.size(), &nvec[0], &svec[0]);
  for (int i = 1; i < svec.size(); i++) ASSERT_TRUE(svec[i] > svec[i-1]);
  if (r2min != 0) { ASSERT_NEAR(svec.front(), 5.36914813e-05, 1e-6); ASSERT_NEAR(svec.back(), 0.685004711, 1e-6); }
  else {            ASSERT_NEAR(svec.front(), 5.23356393e-05, 1e-6); ASSERT_NEAR(svec.back(), 0.682822824, 1e-6); }

  std::vector<float> c0(num_tag, 0.0), c1(num_tag, 0.0), c2(num_tag, 0.0);
  calc.calc_univariate_delta_posterior(trait_index, 0.2, 1.2, 0.1, num_tag, &c0[0], &c1[0], &c2[0]);
  for (int i = 0; i < num_tag; i++) {
    ASSERT_TRUE(c0[i] > 0);
    ASSERT_TRUE(c1[i] != 0);
    ASSERT_TRUE(c2[i] > 0);
    break;
  }

  calc.set_option("fast_cost", 1);
  cost = calc.calc_univariate_cost(trait_index, 0.2, 1.2, 0.1);
  ASSERT_TRUE(std::isfinite(cost));
}

// --gtest_filter=UgmgTest.CalcLikelihood
TEST(UgmgTest, CalcLikelihood) {
  const float r2min = 0.0; 
  const int trait_index = 2; // use second trait for calculations; should work...
  UgmgTest_CalcLikelihood(r2min, trait_index);
}

// --gtest_filter=UgmgTest.CalcLikelihood_with_r2min
TEST(UgmgTest, CalcLikelihood_with_r2min) {
  const float r2min = 0.2;
  const int trait_index = 1;
  UgmgTest_CalcLikelihood(r2min, trait_index);
}

void UgmgTest_CalcLikelihood_testConvolution(float r2min, int trait_index, float pi_val, double costvec[4]) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 10;
  int kmax = 20000; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 1);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", 1);
  calc.set_option("threads", 1);

  calc.set_zvec(trait_index, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait_index, num_tag, &tm.nvec()->at(0));
  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

  const float sig2_beta = 0.1f;
  const float sig2_zeroL = pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;

  calc.set_option("cost_calculator", 0);
  double cost_sampling = calc.calc_univariate_cost(trait_index, pi_val, sig2_zeroA, sig2_beta);
  calc.set_option("cost_calculator", 1);
  double cost_gaussian = calc.calc_univariate_cost(trait_index, pi_val, sig2_zeroA, sig2_beta);
  calc.set_option("cost_calculator", 2);
  double cost_convolve = calc.calc_univariate_cost(trait_index, pi_val, sig2_zeroA, sig2_beta);

  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);
  double cost_unified_gaussian = calc.calc_unified_univariate_cost_gaussian(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL);

  ASSERT_TRUE(std::isfinite(cost_sampling));
  ASSERT_TRUE(std::isfinite(cost_gaussian));
  ASSERT_TRUE(std::isfinite(cost_convolve));
  ASSERT_TRUE(std::isfinite(cost_unified_gaussian));

  // compare that unified gaussian approximation gives the same answer as fast cost function
  // there is a subtle difference here in how do we model inflation arrising form truncated LD structure,
  // therefore with r2min the answer is not precisely the same when r2min!=0.
  if ((r2min==0) || (pi_val==1.0f)) {
    ASSERT_FLOAT_EQ(cost_gaussian, cost_unified_gaussian); 
  }

  std::cout << std::setprecision(9) << cost_sampling << ", " << cost_gaussian << ", " << cost_convolve << ", " << cost_unified_gaussian << std::endl;

  ASSERT_FLOAT_EQ(costvec[0], cost_sampling);
  ASSERT_FLOAT_EQ(costvec[1], cost_gaussian);
  ASSERT_FLOAT_EQ(costvec[2], cost_convolve);
  ASSERT_FLOAT_EQ(costvec[3], cost_unified_gaussian);
}

double calcLikelihoodUnifiedGaussian(float r2min, int trait_index, bool use_complete_tag_indices, float pi_val) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 5;
  int num_tag = 5;  // ideally we should with num_tag < num_snp (i.e. some SNPs having undefined values), but for that it'll be easier to use "init" method
  int N = 100;  // gwas sample size, constant across all variants
  int chr_label = 1;
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", use_complete_tag_indices);
  calc.set_option("threads", 1);

  calc.set_zvec(trait_index, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait_index, num_tag, &tm.nvec()->at(0));
  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(5, &snp_index, &tag_index, &r2);
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure

  const float sig2_beta = 0.1f;
  const float sig2_zeroL = pi_val * sig2_beta;
  const float sig2_zeroA = 1.2f;
  const float sig2_zeroC = 1.0f;

  std::vector<float> pi_vec(num_snp, pi_val);
  std::vector<float> sig2_vec(num_snp, sig2_beta);
  return calc.calc_unified_univariate_cost_gaussian(trait_index, 1, num_snp, &pi_vec[0], &sig2_vec[0], sig2_zeroA, sig2_zeroC, sig2_zeroL);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood
TEST(UgmgTest, CalcConvolveLikelihood) {
  const float r2min = 0.0; 
  const int trait_index = 2; // use second trait for calculations; should work...
  double costvec[4] = {16.0114786, 15.8589964, 15.9299297, 15.8589949};
  UgmgTest_CalcLikelihood_testConvolution(r2min, trait_index, 0.2f, costvec);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihoodInft
TEST(UgmgTest, CalcConvolveLikelihoodInft) {
  const float r2min = 0.0; 
  const int trait_index = 2; // use second trait for calculations; should work...
  double costvec[4] = {1e+100, 20.5189491, 20.5189537, 20.518953};
  UgmgTest_CalcLikelihood_testConvolution(r2min, trait_index, 1.0f, costvec);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood_with_r2min
TEST(UgmgTest, CalcConvolveLikelihood_with_r2min) {
  const float r2min = 0.2;
  const int trait_index = 1;
  double costvec[4] = {16.00285, 15.8589964, 15.9186396, 15.8402427};  
  UgmgTest_CalcLikelihood_testConvolution(r2min, trait_index, 0.2f, costvec);
}

// --gtest_filter=UgmgTest.CalcConvolveLikelihood_with_r2min_inft
TEST(UgmgTest, CalcConvolveLikelihood_with_r2min_inft) {
  const float r2min = 0.2;
  const int trait_index = 1;
  double costvec[4] = {1e+100, 20.5189491, 20.5189472, 20.5189468};  
  UgmgTest_CalcLikelihood_testConvolution(r2min, trait_index, 1.0f, costvec);
}

// --gtest_filter=UgmgTest.CalcUnifiedGaussianLikelihood
TEST(UgmgTest, CalcUnifiedGaussianLikelihood) {
  float v1 = calcLikelihoodUnifiedGaussian(0.0, 1, true, 0.2f);
  float v2 = calcLikelihoodUnifiedGaussian(0.0, 1, false, 0.2f);
  ASSERT_FLOAT_EQ(v1, v2);

  v1 = calcLikelihoodUnifiedGaussian(0.0, 1, true, 1.0f);
  v2 = calcLikelihoodUnifiedGaussian(0.0, 1, false, 1.0f);
  ASSERT_FLOAT_EQ(v1, v2);

  v1 = calcLikelihoodUnifiedGaussian(0.2, 2, true, 0.2f);
  v2 = calcLikelihoodUnifiedGaussian(0.2, 2, false, 0.2f);
  ASSERT_FLOAT_EQ(v1, v2);

  v1 = calcLikelihoodUnifiedGaussian(0.2, 2, true, 1.0f);
  v2 = calcLikelihoodUnifiedGaussian(0.2, 2, false, 1.0f);
  ASSERT_FLOAT_EQ(v1, v2);
}

void BgmgTest_CalcLikelihood_testConvolution(float r2min) {
  // Tests calculation of log likelihood, assuming that all data is already set
  int num_snp = 10;
  int num_tag = 10;
  int kmax = 20000; // #permutations
  int N = 100;  // gwas sample size, constant across all variants
  TestMother tm(num_snp, num_tag, N);
  BgmgCalculator calc;
  calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
  calc.set_option("seed", 0);
  calc.set_option("max_causals", num_snp);
  calc.set_option("kmax", kmax);
  calc.set_option("num_components", 3);
  calc.set_option("cache_tag_r2sum", 1);
  calc.set_option("r2min", r2min);
  calc.set_option("use_complete_tag_indices", 1);
  calc.set_option("threads", 1);

  int trait = 1;
  int chr_label = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  trait = 2; tm.regenerate_zvec();
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);

  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_weights_randprune(20, 0.25);

  float pi_vec[] = { 0.1, 0.2, 0.15 };
  float sig2_beta[] = { 0.5, 0.3 };
  float rho_beta = 0.8;
  float sig2_zero[] = { 1.1, 1.2 };
  float rho_zero = 0.1;

  calc.set_option("cost_calculator", 0);
  double cost_sampling = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  calc.set_option("cost_calculator", 1);
  double cost_gaussian = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  calc.set_option("cost_calculator", 2);
  double cost_convolve = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);

  ASSERT_TRUE(std::isfinite(cost_sampling));
  ASSERT_TRUE(std::isfinite(cost_gaussian));
  ASSERT_TRUE(std::isfinite(cost_convolve));
  std::cout << cost_sampling << ", " << cost_gaussian << ", " << cost_convolve << std::endl;
}

void BgmgTest_CalcLikelihood(float r2min) {
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
  calc.set_option("r2min", r2min);

  int trait = 1;
  int chr_label = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  trait = 2; tm.regenerate_zvec();
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  calc.set_weights(num_tag, &tm.weights()->at(0));

  std::vector<int> snp_index, tag_index;
  std::vector<float> r2;
  tm.make_r2(20, &snp_index, &tag_index, &r2);

  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_weights_randprune(20, 0.25);

  float pi_vec[] = { 0.1, 0.2, 0.15 };
  float sig2_beta[] = { 0.5, 0.3 };
  float rho_beta = 0.8;
  float sig2_zero[] = { 1.1, 1.2 };
  float rho_zero = 0.1;

  double cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  double cost_nocache = calc.calc_bivariate_cost_nocache(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_FLOAT_EQ(cost, cost_nocache);

  // Expect one entry, because calc_bivariate_cost_nocache goes around loglike caching mechanism.
  // If user calls set_option('cache_tag_r2sum', 1), then calc_bivariate_cost will do caching of the log like calculations.
  ASSERT_EQ(calc.get_loglike_cache_size(), 1); int cache_entry = 0;
  calc.get_loglike_cache_bivariate_entry(cache_entry, 3, pi_vec, 2, sig2_beta, &rho_beta, 2, sig2_zero, &rho_zero, &cost);
  ASSERT_FLOAT_EQ(cost, cost_nocache);
  ASSERT_FLOAT_EQ(pi_vec[1], 0.2);
  ASSERT_FLOAT_EQ(sig2_zero[1], 1.2);
  ASSERT_FLOAT_EQ(sig2_beta[0], 0.5);
  ASSERT_FLOAT_EQ(rho_zero, 0.1);
  ASSERT_FLOAT_EQ(rho_beta, 0.8);


  calc.set_option("diag", 0.0);

  calc.set_option("fast_cost", 1);
  cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
  ASSERT_TRUE(std::isfinite(cost));
  ASSERT_EQ(calc.get_loglike_cache_size(), 2);

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

  std::vector<float> c00(num_tag, 0.0), c10(num_tag, 0.0), c01(num_tag, 0.0), c20(num_tag, 0.0), c11(num_tag, 0.0), c02(num_tag, 0.0);
  calc.calc_bivariate_delta_posterior(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, num_tag, &c00[0], &c10[0], &c01[0], &c20[0], &c11[0], &c02[0]);
  for (int i = 0; i < num_tag; i++) {
    ASSERT_TRUE(c00[i] > 0);
    ASSERT_TRUE(c10[i] != 0);
    ASSERT_TRUE(c01[i] != 0);
    ASSERT_TRUE(c20[i] > 0);
    ASSERT_TRUE(c11[i] != 0);
    ASSERT_TRUE(c02[i] > 0);
    break;
  }

  calc.set_option("cache_tag_r2sum", 0);
  calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, zvec_pdf.size(), &zvec1_grid[0], &zvec2_grid[0], &zvec_pdf_nocache[0]);

  for (int i = 0; i < zvec_pdf_nocache.size(); i++)
    ASSERT_FLOAT_EQ(zvec_pdf[i], zvec_pdf_nocache[i]);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood
TEST(BgmgTest, CalcConvolveLikelihood) {
  const float r2min = 0.0; 
  BgmgTest_CalcLikelihood_testConvolution(r2min);
}

// --gtest_filter=BgmgTest.CalcConvolveLikelihood_with_r2min
TEST(BgmgTest, CalcConvolveLikelihood_with_r2min) {
  const float r2min = 0.2;
  BgmgTest_CalcLikelihood_testConvolution(r2min);
}

// bgmg-test.exe --gtest_filter=BgmgTest.CalcLikelihood
TEST(BgmgTest, CalcLikelihood) {
  const float r2min = 0.0;
  BgmgTest_CalcLikelihood(r2min);
}

// bgmg-test.exe --gtest_filter=BgmgTest.CalcLikelihood_with_r2min
TEST(BgmgTest, CalcLikelihood_with_r2min) {
  const float r2min = 0.20;
  BgmgTest_CalcLikelihood(r2min);
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
      int chr_label = 1;
      calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
      calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));
      trait = 2;
      calc.set_zvec(trait, num_tag, &tm2.zvec()->at(0));
      calc.set_nvec(trait, num_tag, &tm2.nvec()->at(0));

      calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
      calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
      calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
      calc.set_ld_r2_csr();  // finalize csr structure
      calc.set_weights_randprune(20, 0.25);

      float pi_vec[] = { 0.1, 0.2, 0.15 };
      float sig2_beta[] = { 0.5, 0.3 };
      float rho_beta = 0.8;
      float sig2_zero[] = { 1.1, 1.2 };
      float rho_zero = 0.1;
      int trait_index = 1;
      double ugmg_cost = calc.calc_univariate_cost(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0]);
      if (num_threads == 1) ugmg_costs[i] = ugmg_cost;
      else ASSERT_FLOAT_EQ(ugmg_costs[i], ugmg_cost);
      ASSERT_FLOAT_EQ(ugmg_cost, calc.calc_univariate_cost(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0]));  // check calc twice => the same cost
      ASSERT_FLOAT_EQ(ugmg_cost, calc.calc_univariate_cost_nocache(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0]));
      
      double bgmg_cost = calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero);
      if (num_threads == 1) bgmg_costs[i] = bgmg_cost;
      else ASSERT_FLOAT_EQ(bgmg_costs[i], bgmg_cost);
      ASSERT_FLOAT_EQ(bgmg_costs[i], calc.calc_bivariate_cost(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero));  // check calc twice => the same cost
      ASSERT_FLOAT_EQ(bgmg_costs[i], calc.calc_bivariate_cost_nocache(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero));

      if ((num_threads == 1) && (i==0)) {
        // test that pdf calculation has the same logic as log likelihood calculation
        std::vector<float> weights(num_tag, 0.0);
        calc.retrieve_weights(num_tag, &weights[0]);

        for (int test_caching = 0; test_caching < 2; test_caching++) {
          calc.set_option("cache_tag_r2sum", test_caching == 0);
          std::vector<float> ugmg_pdf(num_tag, 0.0);
          std::vector<float> bgmg_pdf(num_tag, 0.0);

          for (int tag_pdf_index = 0; tag_pdf_index < num_tag; tag_pdf_index++) {
            std::vector<float> weights2(num_tag, 0.0);
            weights2[tag_pdf_index] = 1.0;
            calc.set_weights(num_tag, &weights2[0]);
            calc.calc_univariate_pdf(trait_index, pi_vec[0], sig2_zero[0], sig2_beta[0], 1, &tm.zvec()->at(tag_pdf_index), &ugmg_pdf[tag_pdf_index]);
            calc.calc_bivariate_pdf(3, pi_vec, 2, sig2_beta, rho_beta, 2, sig2_zero, rho_zero, 1, &tm.zvec()->at(tag_pdf_index), &tm2.zvec()->at(tag_pdf_index), &bgmg_pdf[tag_pdf_index]);
          }
          calc.set_weights(num_tag, &weights[0]); // restore weights back to the original statse

          double ugmg_cost_from_pdf = 0.0, bgmg_cost_from_pdf = 0.0;
          for (int tag_index = 0; tag_index < num_tag; tag_index++) {
            if (weights[tag_index] == 0) continue;
            ugmg_cost_from_pdf += -std::log(static_cast<double>(ugmg_pdf[tag_index])) * weights[tag_index];
            bgmg_cost_from_pdf += -std::log(static_cast<double>(bgmg_pdf[tag_index])) * weights[tag_index];
          }

          ASSERT_FLOAT_EQ(ugmg_cost, ugmg_cost_from_pdf);
          ASSERT_FLOAT_EQ(bgmg_cost, bgmg_cost_from_pdf);
        }
      }
    }

    ASSERT_FLOAT_EQ(ugmg_costs[0], ugmg_costs[2]);
    ASSERT_TRUE(abs(ugmg_costs[0] - ugmg_costs[1]) > 1e-4);
    ASSERT_FLOAT_EQ(bgmg_costs[0], bgmg_costs[2]);
    ASSERT_TRUE(abs(bgmg_costs[0] - bgmg_costs[1]) > 1e-4);
  }
}

void test_tag_r2_caching() {
  int trait_index = 1;
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
  int chr_label = 1;
  calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
  calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

  int sequence_length = 10;
  std::vector<float> num_causal_sequence;
  std::uniform_real_distribution<float> rng(51.0, 99.5);
  for (int i = 0; i < sequence_length; i++) num_causal_sequence.push_back(rng(tm.random_engine()));
  
  calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
  calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
  calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
  calc.set_ld_r2_csr();  // finalize csr structure
  calc.set_weights_randprune(20, 0.25);

  int repeats_count = 100;
  std::vector<double> costs(sequence_length, 0.0);
  for (int j = 0; j < repeats_count; j++) {  // repeat the sequence 100 times and validate that we got the same cost.
    for (int i = 0; i < sequence_length; i++) {
      float pi = static_cast<float>(num_causal_sequence[i]) / static_cast<float>(num_snp);
      double cost = calc.calc_univariate_cost(trait_index, pi, 0.2, 0.15);
      if (j == 0) costs[i] = cost;
      ASSERT_FLOAT_EQ(cost, costs[i]);
    }
  }

  ASSERT_EQ(calc.get_loglike_cache_size(), repeats_count*sequence_length);
  int entry_index = 0;
  for (int j = 0; j < repeats_count; j++) {  // repeat the sequence 100 times and validate that we got the same cost.
    for (int i = 0; i < sequence_length; i++) {
      float pi_vec, sig2_zero, sig2_beta;
      double cost;
      calc.get_loglike_cache_univariate_entry(entry_index++, &pi_vec, &sig2_zero, &sig2_beta, &cost);
      ASSERT_FLOAT_EQ(sig2_zero, 0.2);
      ASSERT_FLOAT_EQ(sig2_beta, 0.15);
      ASSERT_FLOAT_EQ(pi_vec, static_cast<float>(num_causal_sequence[i]) / static_cast<float>(num_snp));
      ASSERT_FLOAT_EQ(cost, costs[i]);
    }
  }
}

// bgmg-test.exe --gtest_filter=Test.tag_r2_caching
TEST(Test, tag_r2_caching) {
  test_tag_r2_caching();
}

// --gtest_filter=Test.performance
TEST(Test, performance) {
  return;
  int trait_index = 1;
  // ideally test with scale = 100. To speedup, use 10 or 1.
  for (int scale = 1; scale <= 100; scale *= 10) {
    SimpleTimer timer_prep(-1);
    std::cout << "scale      : " << scale << "\n";
    int num_snp = 100000 * scale; 
    int num_tag = 10000 * scale;
    int kmax = 10 * scale; // #permutations
    int N = 100000;  // gwas sample size, constant across all variants
    int num_r2 = 10000000 * scale;
    TestMother tm(num_snp, num_tag, N);
    std::vector<int> snp_index, tag_index;
    std::vector<float> r2;
    tm.make_r2(num_r2, &snp_index, &tag_index, &r2);

    BgmgCalculator calc;
    calc.set_tag_indices(num_snp, num_tag, &tm.tag_to_snp()->at(0));
    calc.set_option("max_causals", 0.02 * static_cast<float>(num_snp));
    calc.set_option("kmax", kmax);
    calc.set_option("num_components", 1);
    calc.set_option("cache_tag_r2sum", 1);
    calc.set_option("threads", 32);
    calc.set_seed(123123123);
    int trait = 1;
    int chr_label = 1;
    calc.set_zvec(trait, num_tag, &tm.zvec()->at(0));
    calc.set_nvec(trait, num_tag, &tm.nvec()->at(0));

    calc.set_mafvec(num_snp, &tm.mafvec()->at(0));
    calc.set_chrnumvec(num_snp, &tm.chrnumvec()->at(0));
    calc.set_ld_r2_coo(chr_label, r2.size(), &snp_index[0], &tag_index[0], &r2[0]);
    calc.set_ld_r2_csr();  // finalize csr structure
    calc.set_weights_randprune(20, 0.25);

    std::cout << "preparation: " << timer_prep.elapsed_ms() << " ms\n";
    float pivec[5] = { 0.0001, 0.0003, 0.001, 0.003, 0.01 };
    for (int repeat = 0; repeat < 5; repeat++) {
      SimpleTimer timer(-1);
      double cost_float = calc.calc_univariate_cost_nocache_float(trait_index, pivec[repeat], 0.2, 0.15);
      int time_with_float = timer.elapsed_ms();
      std::cout << "float  cost: " << cost_float << ", time: " << time_with_float << " ms\n";
      ASSERT_TRUE(std::isfinite(cost_float));

      SimpleTimer timer2(-1);
      double cost_double = calc.calc_univariate_cost_nocache_double(trait_index, pivec[repeat], 0.2, 0.15);
      int time_with_double = timer2.elapsed_ms();
      std::cout << "double cost: " << cost_double << ", time: " << time_with_double << " ms\n";
      ASSERT_TRUE(std::isfinite(cost_double));

      std::cout << "cost diff  : " << cost_double - cost_float << "\n";
    }
  }
}

}  // namespace
