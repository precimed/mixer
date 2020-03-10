#pragma once

#include <vector>
#include <string>
#include <map>
#include <algorithm>

#define BED_HEADER_SIZE 3

class BedFileInMemory {
public:
  BedFileInMemory() : num_subjects_(0), num_snps_(0), row_byte_size_(0) {}
  explicit BedFileInMemory(int num_subjects, int num_snps, std::string filename) :
      num_subjects_(num_subjects), num_snps_(num_snps), row_byte_size_((num_subjects + 3) / 4) { read(filename); }

  void read(std::string filename);

  const char* geno(int snp_index) {
    return &buffer_[BED_HEADER_SIZE + snp_index * row_byte_size_];
  }

private:
  const int num_subjects_;
  const int num_snps_;
  const int row_byte_size_;
  std::string buffer_; // remember this has a header (3 bytes), followed by the actual genotypes
};

class BimFile {
public:
  BimFile() {}
  explicit BimFile(std::string filename) { read(filename); }
  explicit BimFile(std::vector<std::string> filenames) { read(filenames); }
  void clear();

  void find_snp_to_index_map();
  int snp_index(const std::string& snp) const;
  int chrposa1a2_index(const std::string& snp) const;
  int chrposa1a2_index(int chri, int bp, const std::string& a1, const std::string& a2) const;

  void read(std::string filename);
  void read(std::vector<std::string> filenames);

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
  std::map<std::string, int> chrposa1a2_to_index_;  // contains all four combinations (chr:bp:a1:a2, chr:bp:a2:a1, and their complements)
};

class FamFile {
public:
  FamFile() {}
  explicit FamFile(std::string filename) { read(filename); }
  void clear();

  void read(std::string filename);

  int size() const { return fid_.size(); }
  const std::vector<std::string>& fid() const { return fid_; }
  const std::vector<std::string>& iid() const { return iid_; }
  const std::vector<std::string>& father_id() const { return father_id_; }
  const std::vector<std::string>& mother_id() const { return mother_id_; }
  const std::vector<int>& sex() const { return sex_; }
  const std::vector<double>& pheno() const { return pheno_; }

private:
  std::vector<std::string> fid_;
  std::vector<std::string> iid_;
  std::vector<std::string> father_id_;
  std::vector<std::string> mother_id_;
  std::vector<int> sex_; // 1 = male, 2 = female, 0 = unknown
  std::vector<double> pheno_;
};

class SnpList {
public:
  SnpList() {}
  explicit SnpList(std::string filename) { read(filename); }
  void read(std::string filename);
  const std::vector<std::string>& snp() const { return snp_; }
  bool contains(const std::string& snp) const;

private:
  std::vector<std::string> snp_;
  std::map<std::string, char> snp_set_;
};

class PlinkLdFile {
public:
  PlinkLdFile(const BimFile& bim, std::string filename);
  void save_as_binary(std::string filename);

private:
  std::vector<int> snpA_index_;
  std::vector<int> snpB_index_;
  std::vector<float> r2_;
};

class FrqFile {
public:
  FrqFile() {}
  FrqFile(const BimFile& bim, std::string filename) { read(bim, filename); }
  FrqFile(const BimFile& bim, std::vector<std::string> filenames) { read(bim, filenames); }
  void read(const BimFile& bim, std::string filename);
  void read(const BimFile& bim, std::vector<std::string> filenames);
  void align_to_reference(const BimFile& bim);
  void clear() { snp_index_.clear(); frq_.clear(); }

  std::vector<float>& frq() { return frq_; }

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
    const std::string& a2reference);

  const std::vector<float>& zscore() { return zscore_; }
  const std::vector<float>& sample_size() { return sample_size_; }

public:
  // Read "N, Z, SNP, A1, A2" 
  // Flip Z scores to align them with the reference
  // SNP     A1      A2      Z       N
  // rs1234567       T       C - 0.154  35217.000
  SumstatFile() {}
  void read(const BimFile& bim, std::string filename);
  void clear() { zscore_.clear(); sample_size_.clear(); }

private:
  std::vector<float> zscore_;
  std::vector<float> sample_size_;
};
