#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <omp.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <future>
#include <map>
#include <limits>

#include <boost/program_options.hpp>
#include <boost/noncopyable.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/filesystem.hpp>
#include <boost/utility.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "bgmg.h"   // VERSION is defined here

namespace po = boost::program_options;
using boost::iostreams::mapped_file_source;

/*
./bgmg_cli.exe --bim-chr H:/GitHub/BGMG/src/11M_bim/chr@.bim.gz --plink-ld H:/work/hapgen_ldmat2_plink/1000Genome_ldmat_p01_SNPwind50k_chr22.ld.gz
./bgmg_cli.exe --bim-chr H:/GitHub/BGMG/src/11M_bim/chr@.bim.gz --frq-chr H:/GitHub/BGMG/src/11M_bim/chr@.frq
*/

class BgmgCpp {
public:
  static void init_log(std::string log_file) {
    bgmg_init_log(log_file.c_str());
  }

  static void log(std::string message) {
    bgmg_log_message(message.c_str());
  }

private:
  int context_id;
};

class Logger : boost::noncopyable {
public:
  Logger() : ss_() {}

  template <typename T>
  std::stringstream& operator<< (const T& rhs) {
    ss_ << rhs;
    return ss_;
  }

  ~Logger() {
    std::cerr << ss_.str() << std::endl;
    BgmgCpp::log(ss_.str());
  }

private:
  std::stringstream ss_;
};

#define LOG Logger()

void log_header(int argc, char *argv[]) {
  std::string header(
    "*********************************************************************\n"
    "* BGMG - Univariate and Bivariate causal mixture models for GWAS     \n"
    "* Version " VERSION "\n"
    "* (C) 2018 Oleksandr Frei et al.,\n"
    "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
    "* GNU General Public License v3\n"
    "*********************************************************************\n");

  LOG << '\n' << header;
  if (argc == 0) return;
  std::stringstream ss;
  ss << "  Call:\n" << argv[0] << " ";
  for (int i = 1; i < argc; i++) {
    if (strlen(argv[i]) == 0) continue;
    if (argv[i][0] == '-') ss << "\\\n\t";
    ss << argv[i] << " ";
  }
  LOG << ss.str();
}

struct BgmgOptions {
  std::string bim;
  std::string bim_chr;
  std::string frq;
  std::string frq_chr;
  std::vector<std::string> bim_files;
  std::vector<std::string> frq_files;
  std::vector<std::string> chr_labels;
  std::string out;
  std::string plink_ld;
  std::string trait1;
};

void describe_bgmg_options(BgmgOptions& s) {
  LOG << "Options in effect (after applying default setting to non-specified parameters):";
  if (!s.bim.empty()) LOG << "\t--bim " << s.bim << " \\";
  if (!s.bim_chr.empty()) LOG << "\t--bim-chr " << s.bim_chr << " \\";
  if (!s.frq.empty()) LOG << "\t--frq " << s.frq << " \\";
  if (!s.frq_chr.empty()) LOG << "\t--frq-chr " << s.frq_chr << " \\";
  if (!s.out.empty()) LOG << "\t--out " << s.out << " \\";
  if (!s.plink_ld.empty()) LOG << "\t--plink-ld " << s.plink_ld << " \\";
  if (!s.trait1.empty()) LOG << "\t--trait1 " << s.trait1 << " \\";
}

class BimFile {
public:
  BimFile() {}
  BimFile(std::string filename) { read(filename); }

  void find_snp_to_index_map() {
    snp_to_index_.clear();
    for (int i = 0; i < size(); i++) {
      if (snp_to_index_.find(snp_[i]) != snp_to_index_.end()) {
        std::stringstream error_str;
        error_str << "Reference contains duplicated variant names (" << snp_[i] << ")";
        throw std::invalid_argument(error_str.str());
      }

      snp_to_index_.insert(std::pair<std::string, int>(snp_[i], i));
    }
  }

  int snp_index(const std::string& snp) const {
    auto iter = snp_to_index_.find(snp);
    if (iter == snp_to_index_.end()) return -1;
    return iter->second;
  }

  void read(std::string filename) {
    //LOG << "Reading " << filename << "...";

    const std::string separators = " \t\n\r";
    std::vector<std::string> tokens;

    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
    int line_no = 0;
    boost::iostreams::filtering_istream in;
    if (boost::algorithm::ends_with(filename, ".gz")) in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    for (std::string str; std::getline(in, str); )
    {
      line_no++;
      int chr_label, bp;
      float gp;
      std::string snp, a1, a2;

      boost::trim_if(str, boost::is_any_of(separators));
      boost::split(tokens, str, boost::is_any_of(separators), boost::token_compress_on);
      try {
        chr_label = stoi(tokens[0]);  // must not use stream or boost::lexical_cast (result in some lock contention)	
        snp = tokens[1];
        gp = stof(tokens[2]);
        bp = stoi(tokens[3]);
        a1 = tokens[4];
        a2 = tokens[5];
      } catch (...) {
        std::stringstream error_str;
        error_str << "Error parsing " << filename << ":" << line_no << " ('" << str << "')";
        throw std::invalid_argument(error_str.str());
      }
      chr_label_.push_back(chr_label);
      snp_.push_back(snp);
      gp_.push_back(gp);
      bp_.push_back(bp);
      a1_.push_back(a1);
      a2_.push_back(a2);
    }

    LOG << "Found " << chr_label_.size() << " variants in " << filename;
  }

  BimFile(std::vector<std::string> filenames) {
    LOG << "Construct reference from " << filenames.size() << " files...";
    std::vector<BimFile> bim_files;
    bim_files.resize(filenames.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < filenames.size(); i++)
      bim_files[i].read(filenames[i]);

    size_t total_size = 0;
    for (int i = 0; i < filenames.size(); i++) total_size += bim_files[i].size();
    chr_label_.reserve(total_size);
    snp_.reserve(total_size);
    gp_.reserve(total_size);
    bp_.reserve(total_size);
    a1_.reserve(total_size);
    a2_.reserve(total_size);

    for (int i = 0; i < filenames.size(); i++) {
      chr_label_.insert(chr_label_.end(), bim_files[i].chr_label_.begin(), bim_files[i].chr_label_.end());
      snp_.insert(snp_.end(), bim_files[i].snp_.begin(), bim_files[i].snp_.end());
      gp_.insert(gp_.end(), bim_files[i].gp_.begin(), bim_files[i].gp_.end());
      bp_.insert(bp_.end(), bim_files[i].bp_.begin(), bim_files[i].bp_.end());
      a1_.insert(a1_.end(), bim_files[i].a1_.begin(), bim_files[i].a1_.end());
      a2_.insert(a2_.end(), bim_files[i].a2_.begin(), bim_files[i].a2_.end());
    }

    LOG << "Found " << chr_label_.size() << " variants in total.";
  }

  int size() const { return chr_label_.size(); }
  const std::vector<int>& chr_label() const { return chr_label_; }
  const std::vector<std::string>& snp() const { return snp_; }
  const std::vector<float>& gp() const { return gp_; }
  const std::vector<int>& bp() const { return bp_; }
  const std::vector<std::string>& a1() const { return a1_; }
  const std::vector<std::string>& a2() const { return a2_; }

private:
  std::vector<int> chr_label_;
  std::vector<std::string> snp_;
  std::vector<float> gp_;
  std::vector<int> bp_;
  std::vector<std::string> a1_;
  std::vector<std::string> a2_;
  std::map<std::string, int> snp_to_index_;
};

class PlinkLdFile {
public:
  PlinkLdFile(const BimFile& bim, std::string filename) {
    const std::string separators = " \t\n\r";
    std::vector<std::string> tokens;

    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
    int line_no = 0;
    int lines_not_match = 0;
    boost::iostreams::filtering_istream in;
    if (boost::algorithm::ends_with(filename, ".gz")) in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    for (std::string str; std::getline(in, str); )
    {
      line_no++;

      if (line_no == 1) continue;  // skip header

      //  CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2
      //    22     16051249   rs62224609     22     16052962  rs376238049     0.774859

      // int chr_a, bp_a, chr_b, bp_b;
      std::string snp_a, snp_b;
      float r2;

      boost::trim_if(str, boost::is_any_of(separators));
      boost::split(tokens, str, boost::is_any_of(separators), boost::token_compress_on);
      try {
        snp_a = tokens[2];
        snp_b = tokens[5];
        r2 = stof(tokens[6]);
      }
      catch (...) {
        std::stringstream error_str;
        error_str << "Error parsing " << filename << ":" << line_no << " ('" << str << "')";
        throw std::invalid_argument(error_str.str());
      }

      int snpA_index = bim.snp_index(snp_a);
      int snpB_index = bim.snp_index(snp_b);
      if (snpA_index < 0 || snpB_index < 0) {
        lines_not_match++;
        continue;
      }

      snpA_index_.push_back(snpA_index);
      snpB_index_.push_back(snpB_index);
      r2_.push_back(r2);

      if (line_no % 100000 == 0) {
        LOG << "Processed " << line_no << " lines";
      }
    }

    LOG << "Parsed " << r2_.size() << " r2 values from " << filename;
    if (lines_not_match) LOG << "[WARNING] " << lines_not_match << " lines ignored because SNP rs# were not found in the reference";
  }
  void save_as_binary(std::string filename) {
    std::ofstream os(filename, std::ofstream::binary);
    if (!os) throw std::runtime_error(::std::runtime_error("can't open" + filename));
    if (sizeof(int) != 4) throw std::runtime_error("sizeof(int) != 4, internal error in BGMG cpp"); // int -> int32_t

    int64_t numel = r2_.size();
    os.write(reinterpret_cast<const char*>(&numel), sizeof(int64_t));

    LOG << "PlinkLdFile::save_as_binary(filename=" << filename << "), writing " << numel << " elements...";
    os.write(reinterpret_cast<char*>(&snpA_index_[0]), numel * sizeof(int));
    os.write(reinterpret_cast<char*>(&snpB_index_[0]), numel * sizeof(int));
    os.write(reinterpret_cast<char*>(&r2_[0]), numel * sizeof(float));
    os.close();
  }
private:
  std::vector<int> snpA_index_;
  std::vector<int> snpB_index_;
  std::vector<float> r2_;
};

class FrqFile {
public:
  FrqFile() {
  }
  FrqFile(const BimFile& bim, std::string filename) {
    read(bim, filename);
  }

  void read(const BimFile& bim, std::string filename) {
    const std::string separators = " \t\n\r";
    std::vector<std::string> tokens;

    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
    int line_no = 0;
    int lines_not_match = 0;
    boost::iostreams::filtering_istream in;
    if (boost::algorithm::ends_with(filename, ".gz")) in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    for (std::string str; std::getline(in, str); )
    {
      line_no++;

      if (line_no == 1) continue;  // skip header

                                   //  CHR           SNP   A1   A2          MAF  NCHROBS
                                   // 1   rs376342519 CCGCCGTTGCAAAGGCGCGCCG    C     0.005964     1006

      std::string snp;
      float frq;

      boost::trim_if(str, boost::is_any_of(separators));
      boost::split(tokens, str, boost::is_any_of(separators), boost::token_compress_on);
      try {
        snp = snp = tokens[1];
        frq = stof(tokens[4]);
      }
      catch (...) {
        std::stringstream error_str;
        error_str << "Error parsing " << filename << ":" << line_no << " ('" << str << "')";
        throw std::invalid_argument(error_str.str());
      }

      int snp_index = bim.snp_index(snp);
      if (snp_index < 0) {
        lines_not_match++;
        continue;
      }

      snp_index_.push_back(snp_index);
      frq_.push_back(frq);
    }

    LOG << "Parsed " << frq_.size() << " frq values from " << filename;
    if (lines_not_match) LOG << "[WARNING] " << lines_not_match << " lines ignored because SNP rs# were not found in the reference";
  }

  FrqFile(const BimFile& bim, std::vector<std::string> filenames) {
    LOG << "Construct frq from " << filenames.size() << " files...";
    std::vector<FrqFile> frq_files;
    frq_files.resize(filenames.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < filenames.size(); i++)
      frq_files[i].read(bim, filenames[i]);

    size_t total_size = 0;
    for (int i = 0; i < filenames.size(); i++) total_size += frq_files[i].frq_.size();
    snp_index_.reserve(total_size);
    frq_.reserve(total_size);

    for (int i = 0; i < filenames.size(); i++) {
      snp_index_.insert(snp_index_.end(), frq_files[i].snp_index_.begin(), frq_files[i].snp_index_.end());
      frq_.insert(frq_.end(), frq_files[i].frq_.begin(), frq_files[i].frq_.end());
    }

    LOG << "Found " << snp_index_.size() << " frq values in total.";
  }

  void align_to_reference(const BimFile& bim) {
    std::vector<float> new_frq(bim.size(), std::numeric_limits<float>::infinity());
    for (int i = 0; i < frq_.size(); i++) {
      new_frq[snp_index_[i]] = frq_[i];
    }
    frq_.swap(new_frq);
    snp_index_.swap(std::vector<int>());

    for (int i = 0; i < frq_.size(); i++) {
      if (!std::isfinite(frq_[i])) {
        std::stringstream error_str;
        error_str << "Variant " << bim.snp()[i] << " is not present in .frq files";
        throw std::runtime_error(error_str.str());
      }
    }
  }

private:
  std::vector<int> snp_index_;
  std::vector<float> frq_;
};

class SumstatFile {
 public:
  enum FLIP_STATUS {
    FLIP_STATUS_ALIGNED,
    FLIP_STATUS_FLIPPED,
    FLIP_STATUS_AMBIGUOUS,
    FLIP_STATUS_MISMATCH,
  };

  // Detect flip status
  static FLIP_STATUS flip_strand(
    const std::string& a1sumstat,
    const std::string& a2sumstat,
    const std::string& a1reference,
    const std::string& a2reference) {

    std::string a1 = boost::to_lower_copy(a1sumstat);
    std::string a2 = boost::to_lower_copy(a2sumstat);
    std::string a1ref = boost::to_lower_copy(a1reference);
    std::string a2ref = boost::to_lower_copy(a2reference);

    // validate that all charactars are A, T, C, G - nothing else
    auto check_atcg = [](std::string val) {
      for (int i = 0; i < val.size(); i++) { if (val[i] != 'a' && val[i] != 't' && val[i] != 'c' && val[i] != 'g') return false; }
      return true;
    };

    auto atcg_complement = [](std::string val) {
      std::string retval;
      for (int i = 0; i < val.size(); i++) {
        if (val[i] == 'a') retval += 't';
        if (val[i] == 't') retval += 't';
        if (val[i] == 'c') retval += 'g';
        if (val[i] == 'g') retval += 'c';
        return retval;
      }
    };

    if (!check_atcg(a1 + a2 + a1ref + a2ref)) return FLIP_STATUS_MISMATCH;
    if (atcg_complement(a1) == a2) return FLIP_STATUS_AMBIGUOUS;

    if (a1 == a1ref && a2 == a2ref) return FLIP_STATUS_ALIGNED;
    if (a1 == a2ref && a2 == a1ref) return FLIP_STATUS_FLIPPED;
    if (atcg_complement(a1) == a1ref && atcg_complement(a2) == a2ref) return FLIP_STATUS_FLIPPED;
    if (atcg_complement(a1) == a2ref && atcg_complement(a2) == a1ref) return FLIP_STATUS_ALIGNED;

    return FLIP_STATUS_MISMATCH;
  }

 public:
  // Read "N, Z, SNP, A1, A2" 
  // Flip Z scores to align them with the reference
  // SNP     A1      A2      Z       N
  // rs1234567       T       C - 0.154  35217.000
  SumstatFile(const BimFile& bim, std::string filename) {
    std::vector<int> snp_index_;
    const std::string separators = " \t\n\r";
    std::vector<std::string> tokens;

    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);

    // gather statistics
    int line_no = 0;
    int lines_incomplete = 0;
    int snps_dont_match_reference = 0;
    int snp_ambiguous = 0;
    int mismatch_alleles = 0;
    int flipped_alleles = 0;
    int duplicates_ignored = 0;
    
    boost::iostreams::filtering_istream in;
    if (boost::algorithm::ends_with(filename, ".gz")) in.push(boost::iostreams::gzip_decompressor());
    in.push(file);

    int snp_col = -1, a1_col = -1, a2_col = -1, z_col = -1, n_col = -1;
    int num_cols_required;
    for (std::string str; std::getline(in, str); )
    {
      line_no++;
      boost::trim_if(str, boost::is_any_of(separators));
      boost::to_lower(str);
      boost::split(tokens, str, boost::is_any_of(separators), boost::token_compress_on);

      if (line_no == 1) { // process header
        for (int coli = 0; coli < tokens.size(); coli++) {
          if (tokens[coli] == std::string("snp")) snp_col = coli;
          if (tokens[coli] == std::string("a1")) a1_col = coli;
          if (tokens[coli] == std::string("a2")) a2_col = coli;
          if (tokens[coli] == std::string("n")) n_col = coli;
          if (tokens[coli] == std::string("z")) z_col = coli;
        }

        if (snp_col == -1 || a1_col == -1 || a2_col == -1 || z_col == -1 || n_col == -1) {
          std::stringstream error_str;
          error_str << "Error parsing " << filename << ", unexpected header: " << str;
          throw std::invalid_argument(error_str.str());
        }

        num_cols_required = 1 + std::max({ snp_col, a1_col, a2_col, n_col, z_col });

        continue;  // finish processing header
      }

      if (tokens.size() < num_cols_required) {
        lines_incomplete++;
        continue;
      }

      std::string snp, a1, a2;
      float sample_size, zscore;

      try {
        snp = tokens[snp_col];
        a1 = tokens[a1_col];
        a2 = tokens[a2_col];
        zscore = stof(tokens[z_col]);
        sample_size = stof(tokens[n_col]);
      }
      catch (...) {
        std::stringstream error_str;
        error_str << "Error parsing " << filename << ":" << line_no << " ('" << str << "')";
        throw std::invalid_argument(error_str.str());
      }

      int snp_index = bim.snp_index(snp);
      if (snp_index < 0) {
        snps_dont_match_reference++;
        continue;
      }

      FLIP_STATUS flip_status = flip_strand(a1, a2, bim.a1()[snp_index], bim.a2()[snp_index]);
      if (flip_status == FLIP_STATUS_AMBIGUOUS) {
        snp_ambiguous++;
        continue;
      }
      if (flip_status == FLIP_STATUS_MISMATCH) {
        mismatch_alleles++;
        continue;
      }

      if (flip_status == FLIP_STATUS_FLIPPED) {
        flipped_alleles++;
        zscore *= -1.0f;
      }

      zscore_.push_back(zscore);
      sample_size_.push_back(sample_size);
      snp_index_.push_back(snp_index);
    }

    // align to reference
    {
      std::vector<float> new_zscore(bim.size(), std::numeric_limits<float>::infinity());
      std::vector<float> new_sample_size(bim.size(), std::numeric_limits<float>::infinity());
      for (int i = 0; i < snp_index_.size(); i++) {
        if (std::isfinite(new_zscore[snp_index_[i]]) || std::isfinite(new_sample_size[snp_index_[i]])) {
          duplicates_ignored++;
          continue;
        }

        new_zscore[snp_index_[i]] = zscore_[i];
        new_sample_size[snp_index_[i]] = sample_size_[i];
      }
      zscore_.swap(new_zscore);
      sample_size_.swap(new_sample_size);

      int num_def = 0;
      for (int i = 0; i < zscore_.size(); i++) {
        if (std::isfinite(zscore_[i]) && std::isfinite(sample_size_[i])) num_def++;
      }

      LOG << "Found " << num_def << " variants with well-defined Z and N in " << filename << ". Other statistics: ";
      LOG << "\t" << line_no << " lines found (including header)";
      if (lines_incomplete > 0) LOG << "\t" << lines_incomplete << " lines were ignored as there are incomplete (too few values).";
      if (duplicates_ignored > 0) LOG << "\t" << duplicates_ignored << " lines were ignored as they contain duplicated RS#.";
      if (snps_dont_match_reference > 0) LOG << "\t" << snps_dont_match_reference << " lines were ignored as RS# does not match reference file.";
      if (mismatch_alleles > 0) LOG << "\t" << mismatch_alleles << " variants were ignored as they had A1/A2 alleles that do not match reference.";
      if (snp_ambiguous > 0) LOG << "\t" << snp_ambiguous << " variants were ignored as they are strand-ambiguous.";
      if (flipped_alleles > 0) LOG << "\t" << flipped_alleles << " variants had flipped A1/A2 alleles; sign of z-score was flipped.";
    }
  }

 private:
   std::vector<float> zscore_;
   std::vector<float> sample_size_;
};

void fix_and_validate(BgmgOptions& bgmg_options, po::variables_map& vm) {
  // initialize chr_labels
  if (bgmg_options.chr_labels.empty()) {
    for (int i = 1; i <= 22; i++)
      bgmg_options.chr_labels.push_back(boost::lexical_cast<std::string>(i));
  }

  // Validate --bim / --bim-chr option
  if (bgmg_options.bim.empty() && bgmg_options.bim_chr.empty())
    throw std::invalid_argument(std::string("ERROR: Either --bim or --bim-chr must be specified"));

  if (!bgmg_options.bim_chr.empty()) {
    if (!boost::contains(bgmg_options.bim_chr, "@"))
      throw std::invalid_argument(std::string("ERROR: --bim-chr must have @ to indicate location of chr label"));

    for (auto chrlabel : bgmg_options.chr_labels) {
      bgmg_options.bim_files.push_back(bgmg_options.bim_chr);
      boost::replace_all(bgmg_options.bim_files.back(), "@", chrlabel);
    }
  }
  else {
    bgmg_options.bim_files.push_back(bgmg_options.bim);
  }

  for (auto& bim_file: bgmg_options.bim_files)
  if (!boost::filesystem::exists(bim_file)) {
    std::stringstream ss; ss << "ERROR: input file " << bim_file << " does not exist";
    throw std::runtime_error(ss.str());
  }

  // Validate --out option
  if (bgmg_options.out.empty())
    throw std::invalid_argument(std::string("ERROR: --out option must be specified"));

  // Validate --plink-ld option, and stop further validation if plink_ld is enabled.
  if (!bgmg_options.plink_ld.empty()) {
    if (!boost::filesystem::exists(bgmg_options.plink_ld)) {
      std::stringstream ss; ss << "ERROR: input file " << bgmg_options.plink_ld << " does not exist";
      throw std::runtime_error(ss.str());
    }

    return;
  }

  // Validate --frq / --frq-chr option
  if (bgmg_options.frq.empty() && bgmg_options.frq_chr.empty())
    throw std::invalid_argument(std::string("ERROR: Either --frq or --frq-chr must be specified"));

  if (!bgmg_options.frq_chr.empty()) {
    if (!boost::contains(bgmg_options.frq_chr, "@"))
      throw std::invalid_argument(std::string("ERROR: --frq-chr must have @ to indicate location of chr label"));

    for (auto chrlabel : bgmg_options.chr_labels) {
      bgmg_options.frq_files.push_back(bgmg_options.frq_chr);
      boost::replace_all(bgmg_options.frq_files.back(), "@", chrlabel);
    }
  }
  else {
    bgmg_options.frq_files.push_back(bgmg_options.frq);
  }

  for (auto& frq_file : bgmg_options.frq_files) {
    if (!boost::filesystem::exists(frq_file)) {
      std::stringstream ss; ss << "ERROR: input file " << frq_file << " does not exist";
      throw std::runtime_error(ss.str());
    }
  }

  // Validate trait1 option
  if (bgmg_options.trait1.empty() || !boost::filesystem::exists(bgmg_options.trait1))
    throw std::invalid_argument(std::string("ERROR: Either --trait1 file does not exist: " + bgmg_options.trait1));
}

int main(int argc, char *argv[]) {
  try {
    BgmgOptions bgmg_options;
    po::options_description po_options("BGMG " VERSION " - Univariate and Bivariate causal mixture models for GWAS");
    po_options.add_options()
      ("help,h", "produce this help message")
      ("bim", po::value(&bgmg_options.bim), "Path to .bim file that defines the reference set of SNPs")
      ("bim-chr", po::value(&bgmg_options.bim_chr), "Path to .bim files, split per chromosome. Use @ symbol to indicate location of chromosome label.")
      ("frq", po::value(&bgmg_options.frq), "Path to .frq file that defines the minor allele frequency for the reference set of SNPs")
      ("frq-chr", po::value(&bgmg_options.frq_chr), "Path to .frq files, split per chromosome. Use @ symbol to indicate location of chromosome label.")
      ("plink-ld", po::value(&bgmg_options.plink_ld), "Path to plink .ld.gz file to convert into BGMG binary format.")
      ("chr-labels", po::value(&bgmg_options.chr_labels)->multitoken(), "Set of chromosome labels. Defaults to '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'")
      ("out", po::value(&bgmg_options.out)->default_value("bgmg"),
        "prefix of the output files; "
        "See README.md file for detailed description of file formats.")
      ("trait1", po::value(&bgmg_options.trait1), "Path to .sumstats.gz file for the trait to analyze")
    ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(po_options)
      .run(), vm);
    notify(vm);

    bool show_help = (vm.count("help") > 0);
    if (show_help) {
      std::cerr << po_options;
      exit(EXIT_FAILURE);
    }

    BgmgCpp::init_log(bgmg_options.out + ".bgmglib.log");
    log_header(argc, argv);

    try {
      auto analysis_started = boost::posix_time::second_clock::local_time();
      LOG << "Analysis started: " << analysis_started;
      fix_and_validate(bgmg_options, vm);
      describe_bgmg_options(bgmg_options);

      BimFile bim_file(bgmg_options.bim_files);
      bim_file.find_snp_to_index_map();

      if (!bgmg_options.plink_ld.empty()) {
        PlinkLdFile plink_ld_file(bim_file, bgmg_options.plink_ld);
        plink_ld_file.save_as_binary(bgmg_options.out + ".ld.bin");
      } else {
        FrqFile frq_file(bim_file, bgmg_options.frq_files);
        frq_file.align_to_reference(bim_file);

        SumstatFile trait1_file(bim_file, bgmg_options.trait1);

        // ready to extract zvec, nvec, mafvec, chrnumvec, posvec (indexed as the reference, not as tag SNPs)
      }

      auto analysis_finished = boost::posix_time::second_clock::local_time();
      LOG << "Analysis finished: " << analysis_finished;
      LOG << "Elapsed time: " << analysis_finished - analysis_started;
    }
    catch (std::exception& e) {
      LOG << "ERROR: " << e.what();
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << "Exception  : " << e.what() << "\n";
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown error occurred.";
    return EXIT_FAILURE;
  }
}
