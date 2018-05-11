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

#include "bgmg_calculator.h"

int64_t BgmgCalculator::set_zvec(int trait, int length, double* values) {
  return 0;
}

int64_t BgmgCalculator::set_nvec(int trait, int length, double* values) {
  return 0;
}

int64_t BgmgCalculator::set_hvec(int length, double* values) {
  return 0;
}

int64_t BgmgCalculator::set_w_ld(int length, double* values) {
  return 0;
}

int64_t BgmgCalculator::set_ref_ld_sum_r2(int length, int r2bins, double* values) {
  return 0;
}

int64_t BgmgCalculator::set_ref_ld_sum_r4(int length, int r2bins, double* values) {
  return 0;
}

int64_t BgmgCalculator::set_option(char* option, int value) {
  return 0;
}

double BgmgCalculator::calc_univariate_cost(double pi_vec, double sig2_zero, double sig2_beta) {
  return 0.0;
}

double BgmgCalculator::calc_bivariate_cost(int num_components, int num_traits, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero) {
  return 0.0;
}
