#include "bgmg_parse.h"

#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "bgmg_log.h"
#include "zstr.hpp"


void BimFile::find_snp_to_index_map() {
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

int BimFile::snp_index(const std::string& snp) const {
  auto iter = snp_to_index_.find(snp);
  if (iter == snp_to_index_.end()) return -1;
  return iter->second;
}

std::shared_ptr<std::istream> open_file(std::string filename) {
  if (boost::algorithm::ends_with(filename, ".gz")) {
    return std::make_shared<zstr::ifstream>(filename, std::ios_base::in | std::ios_base::binary);
  }
  else {
    return std::make_shared<std::ifstream>(filename, std::ios_base::in | std::ios_base::binary);
  }
}

void BimFile::read(std::string filename) {
  //LOG << "Reading " << filename << "...";

  const std::string separators = " \t\n\r";
  std::vector<std::string> tokens;

  std::shared_ptr<std::istream> in_ptr = open_file(filename);
  std::istream& in = *in_ptr;

  int line_no = 0;
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
    }
    catch (...) {
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

  LOG << " Found " << chr_label_.size() << " variants in " << filename;
}

void BimFile::read(std::vector<std::string> filenames) {
  LOG << " Construct reference from " << filenames.size() << " files...";
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

  LOG << " Found " << chr_label_.size() << " variants in total.";
}

PlinkLdFile::PlinkLdFile(const BimFile& bim, std::string filename) {
  const std::string separators = " \t\n\r";
  std::vector<std::string> tokens;

  LOG << " Reading " << filename << "...";

  std::shared_ptr<std::istream> in_ptr = open_file(filename);
  std::istream& in = *in_ptr;
  int line_no = 0;
  int lines_not_match = 0;
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
      LOG << " Processed " << line_no << " lines";
    }
  }

  LOG << " Parsed " << r2_.size() << " r2 values from " << filename;
  if (lines_not_match) LOG << " [WARNING] " << lines_not_match << " lines ignored because SNP rs# were not found in the reference";
}

void PlinkLdFile::save_as_binary(std::string filename) {
  std::ofstream os(filename, std::ofstream::binary);
  if (!os) throw std::runtime_error(::std::runtime_error("can't open" + filename));
  if (sizeof(int) != 4) throw std::runtime_error("sizeof(int) != 4, internal error in BGMG cpp"); // int -> int32_t

  int64_t numel = r2_.size();
  os.write(reinterpret_cast<const char*>(&numel), sizeof(int64_t));

  LOG << " PlinkLdFile::save_as_binary(filename=" << filename << "), writing " << numel << " elements...";
  os.write(reinterpret_cast<char*>(&snpA_index_[0]), numel * sizeof(int));
  os.write(reinterpret_cast<char*>(&snpB_index_[0]), numel * sizeof(int));
  os.write(reinterpret_cast<char*>(&r2_[0]), numel * sizeof(float));
  os.close();
}


void FrqFile::read(const BimFile& bim, std::string filename) {
  const std::string separators = " \t\n\r";
  std::vector<std::string> tokens;

  std::shared_ptr<std::istream> in_ptr = open_file(filename);
  std::istream& in = *in_ptr;
  int line_no = 0;
  int lines_not_match = 0;
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
      snp = tokens[1];
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

  LOG << " Parsed " << frq_.size() << " frq values from " << filename;
  if (lines_not_match) LOG << " [WARNING] " << lines_not_match << " lines ignored because SNP rs# were not found in the reference";
}

void FrqFile::read(const BimFile& bim, std::vector<std::string> filenames) {
  if (filenames.empty()) return;

  LOG << " Construct frq from " << filenames.size() << " files...";
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

  LOG << " Found " << snp_index_.size() << " frq values in total.";
}

void FrqFile::align_to_reference(const BimFile& bim) {
  std::vector<float> new_frq(bim.size(), std::numeric_limits<float>::infinity());
  for (int i = 0; i < frq_.size(); i++) {
    new_frq[snp_index_[i]] = frq_[i];
  }
  frq_.swap(new_frq);
  snp_index_.clear();

  for (int i = 0; i < frq_.size(); i++) {
    if (!std::isfinite(frq_[i])) {
      std::stringstream error_str;
      error_str << "Variant " << bim.snp()[i] << " is not present in .frq files";
      throw std::runtime_error(error_str.str());
    }
  }
}

// Detect flip status
SumstatFile::FLIP_STATUS SumstatFile::flip_strand(
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
    }
    return retval;
  };

  if (!check_atcg(a1 + a2 + a1ref + a2ref)) return FLIP_STATUS_MISMATCH;
  if (atcg_complement(a1) == a2) return FLIP_STATUS_AMBIGUOUS;

  if (a1 == a1ref && a2 == a2ref) return FLIP_STATUS_ALIGNED;
  if (a1 == a2ref && a2 == a1ref) return FLIP_STATUS_FLIPPED;
  if (atcg_complement(a1) == a1ref && atcg_complement(a2) == a2ref) return FLIP_STATUS_FLIPPED;
  if (atcg_complement(a1) == a2ref && atcg_complement(a2) == a1ref) return FLIP_STATUS_ALIGNED;

  return FLIP_STATUS_MISMATCH;
}

// Read "N, Z, SNP, A1, A2" 
// Flip Z scores to align them with the reference
// SNP     A1      A2      Z       N
// rs1234567       T       C - 0.154  35217.000
void SumstatFile::read(const BimFile& bim, std::string filename) {
  std::vector<int> snp_index_;
  const std::string separators = " \t\n\r";
  std::vector<std::string> tokens;

  std::shared_ptr<std::istream> in_ptr = open_file(filename);
  std::istream& in = *in_ptr;

  // gather statistics
  int line_no = 0;
  int lines_incomplete = 0;
  int snps_dont_match_reference = 0;
  int snp_ambiguous = 0;
  int mismatch_alleles = 0;
  int flipped_alleles = 0;
  int duplicates_ignored = 0;
  int undefined_z_or_n = 0;

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

    if ((tokens[z_col] == "na") || (tokens[n_col] == "na")) {
      undefined_z_or_n++;
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

    LOG << " Found " << num_def << " variants with well-defined Z and N in " << filename << ". Other statistics: ";
    LOG << " \t" << line_no << " lines found (including header)";
    if (lines_incomplete > 0) LOG << " \t" << lines_incomplete << " lines were ignored as there are incomplete (too few values).";
    if (undefined_z_or_n > 0) LOG << " \t" << undefined_z_or_n << " lines were ignored as Z or N value was not define ('na' or similar).";
    if (duplicates_ignored > 0) LOG << " \t" << duplicates_ignored << " lines were ignored as they contain duplicated RS#.";
    if (snps_dont_match_reference > 0) LOG << " \t" << snps_dont_match_reference << " lines were ignored as RS# does not match reference file.";
    if (mismatch_alleles > 0) LOG << " \t" << mismatch_alleles << " variants were ignored as they had A1/A2 alleles that do not match reference.";
    if (snp_ambiguous > 0) LOG << " \t" << snp_ambiguous << " variants were ignored as they are strand-ambiguous.";
    if (flipped_alleles > 0) LOG << " \t" << flipped_alleles << " variants had flipped A1/A2 alleles; sign of z-score was flipped.";
  }
}

void BimFile::clear() {
  chr_label_.clear();
  snp_.clear();
  gp_.clear();
  bp_.clear();
  a1_.clear();
  a2_.clear();
  snp_to_index_.clear();
}

void SnpList::read(std::string filename) {
  const std::string separators = " \t\n\r";
  std::shared_ptr<std::istream> in_ptr = open_file(filename);
  std::istream& in = *in_ptr;
  std::stringstream buffer;
  buffer << in.rdbuf();
  std::string str = buffer.str();
  boost::trim_if(str, boost::is_any_of(separators));
  boost::split(snp_, str, boost::is_any_of(separators), boost::token_compress_on);

  snp_set_.clear();
  for (int i = 0; i < snp_.size(); i++) {
    if (!snp_[i].empty()) {
      snp_set_.insert(std::pair<std::string, char>(boost::to_lower_copy(snp_[i]), 1));
    }
  }
}

bool SnpList::contains(const std::string& snp) const { 
  return snp_set_.find(boost::to_lower_copy(snp)) != snp_set_.end();
}

void FamFile::clear() {
  fid_.clear();
  iid_.clear();
  father_id_.clear();
  mother_id_.clear();
  sex_.clear();
  pheno_.clear();
}

void FamFile::read(std::string filename) {
//LOG << "Reading " << filename << "...";

  const std::string separators = " \t\n\r";
  std::vector<std::string> tokens;

  std::shared_ptr<std::istream> in_ptr = open_file(filename);
  std::istream& in = *in_ptr;

  int line_no = 0;
  for (std::string str; std::getline(in, str); )
  {
    line_no++;
    std::string fid, iid, father_id, mother_id;
    int sex;
    double pheno;

    boost::trim_if(str, boost::is_any_of(separators));
    boost::split(tokens, str, boost::is_any_of(separators), boost::token_compress_on);
    try {
      fid = tokens[0];
      iid = tokens[1];
      father_id = tokens[2];
      mother_id = tokens[3];
      sex = stoi(tokens[4]);  // must not use stream or boost::lexical_cast (result in some lock contention)
      pheno = stod(tokens[5]);
    }
    catch (...) {
      std::stringstream error_str;
      error_str << "Error parsing " << filename << ":" << line_no << " ('" << str << "')";
      throw std::invalid_argument(error_str.str());
    }
    fid_.push_back(fid);
    iid_.push_back(iid);
    father_id_.push_back(father_id);
    mother_id_.push_back(mother_id);
    sex_.push_back(sex);
    pheno_.push_back(pheno);
  }

  LOG << " Found " << fid_.size() << " individuals in " << filename;
}

void BedFileInMemory::read(std::string filename) {
  LOG << " Reading " << num_subjects_ << " subjects, " << num_snps_ << " variants from " << filename;
  std::ifstream bedfile_stream(filename, std::ios::in | std::ios::binary);
  std::istreambuf_iterator<char> eos;
  std::string buffer(std::istreambuf_iterator<char>(bedfile_stream), eos);
  buffer_.swap(buffer);
  if (buffer_.size() != row_byte_size_ * num_snps_ + BED_HEADER_SIZE) throw std::invalid_argument("plink .bed file has a wrong size");
  LOG << " Finish reading " << filename;
}