#pragma once

#include <vector>
#include <string>
#include <map>

class BimFile {
public:
  BimFile() {}
  explicit BimFile(std::string filename) { read(filename); }
  explicit BimFile(std::vector<std::string> filenames) { read(filenames); }
  void clear();

  void find_snp_to_index_map();
  int snp_index(const std::string& snp) const;
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
