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

#define VERSION "v0.9.0"

#include <stdexcept>
#include "boost/exception/diagnostic_information.hpp"
#include "boost/exception/get_error_info.hpp"

#include <bgmg.h>
#include "bgmg_calculator.h"

static std::string last_error_;
static const char* last_error() { return last_error_.c_str(); }
static void set_last_error(const std::string& error) { last_error_.assign(error); }
const char* bgmg_get_last_error() { return last_error(); }

#define CATCH_EXCEPTIONS                                                       \
catch (const std::runtime_error& e) {                                          \
  set_last_error("runtime_error:  " + std::string(e.what()));                  \
  return -1;                                                                   \
} catch (...) {                                                                \
  set_last_error(boost::current_exception_diagnostic_information());           \
  return -1;                                                                   \
}

int64_t bgmg_set_zvec(int context_id, int trait, int length, double* values) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->set_zvec(trait, length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_nvec(int context_id, int trait, int length, double* values) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->set_nvec(trait, length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_hvec(int context_id, int length, double* values) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->set_hvec(length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_weights(int context_id, int length, double* values) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->set_weights(length, values);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_ref_ld(int context_id, int length, int r2bins, double* sum_r2, double* sum_r4) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->set_ref_ld(length, r2bins, sum_r2, sum_r4);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_set_option(int context_id, char* option, int value) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->set_option(option, value);
  } CATCH_EXCEPTIONS;
}

double bgmg_calc_univariate_cost(int context_id, double pi_vec, double sig2_zero, double sig2_beta) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_univariate_cost(pi_vec, sig2_zero, sig2_beta);
  } CATCH_EXCEPTIONS;
}

double bgmg_calc_bivariate_cost(int context_id, int num_components, double* pi_vec, double* sig2_beta, double* rho_beta, double* sig2_zero, double rho_zero) {
  try {
    return BgmgCalculatorManager::singleton().Get(context_id)->calc_bivariate_cost(num_components, pi_vec, sig2_beta, rho_beta, sig2_zero, rho_zero);
  } CATCH_EXCEPTIONS;
}

int64_t bgmg_dispose(int context_id) {
  try {
    BgmgCalculatorManager::singleton().Erase(context_id);
    return 0;
  } CATCH_EXCEPTIONS;
}
