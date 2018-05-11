/*
  bgmg - tool to calculate log likelihood of BGMG and UGMG mixture models
  Copyright (C) 2018 Oleksandr Frei 

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <stdint.h>

#include <unordered_map>
#include <memory>

#include "boost/utility.hpp"

// Singleton class to manage a collection of objects, identifiable with some integer ID.
template<class Type>
class TemplateManager : boost::noncopyable {
 public:
  static TemplateManager<Type>& singleton() {
    static TemplateManager<Type> manager;
    return manager;
  }

  std::shared_ptr<Type> Get(int id) {
    if (map_.find(id) == map_.end()) {
      map_[id] = std::make_shared<Type>();
    }
    return map_[id];
  }

  void Erase(int id) {
    map_.erase(id);
  }

 private:
  TemplateManager() { }  // Singleton (make constructor private)
  std::unordered_map<int, std::shared_ptr<Type>> map_;
};


class BgmgCalculator {
 public:
  BgmgCalculator() : num_snps_(-1), r2bins_(-1), k_max_(100) {}
  int64_t set_zvec(int trait, int length, double* values);
  int64_t set_nvec(int trait, int length, double* values);
  int64_t set_hvec(int length, double* values);
  int64_t set_weights(int length, double* values);
  int64_t set_ref_ld(int length, int r2bins, double* sum_r2, double* sum_r4);
  int64_t set_option(char* option, int value);
  double calc_univariate_cost(float pi_vec, float sig2_zero, float sig2_beta);
  double calc_bivariate_cost(int num_components, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero);
 private:
  int num_snps_;
  int r2bins_;
  std::vector<float> zvec1_;
  std::vector<float> nvec1_;
  std::vector<float> hvec_;
  std::vector<float> w_;
  std::vector<float> r2_;       // matrix of size num_snps * r2bins, containing effective r2 in each bin
  std::vector<int>   r2_hist_;  // matrix of size num_snps * r2bins, containing an integer number of SNPs in that bin
  int k_max_;
  void set_num_snps(int length);
};

typedef TemplateManager<BgmgCalculator> BgmgCalculatorManager;
