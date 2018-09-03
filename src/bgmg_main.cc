#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <omp.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <future>

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

// Copied from https://github.com/bigartm/bigartm/blob/master/src/artm/utility/ifstream_or_cin.h
class ifstream_or_cin {
public:
  explicit ifstream_or_cin(const std::string& filename) {
    if (filename == "-")  // read from std::cin
      return;

    if (!boost::filesystem::exists(filename))
      throw std::runtime_error("File " + filename + " does not exist.");

    if (boost::filesystem::exists(filename) && !boost::filesystem::is_regular_file(filename))
      throw std::runtime_error("File " + filename + " is not regular (probably it's a directory).");

    file_.open(filename);
  }

  std::istream& get_stream() { return file_.is_open() ? file_ : std::cin; }

  size_t size() {
    if (!file_.is_open()) {
      return 0;
    }

    size_t pos = file_.tellg();
    file_.seekg(0, std::ios::end);
    size_t size = file_.tellg();
    file_.seekg(pos, std::ios::beg);
    return size;
  }

private:
  boost::iostreams::stream<mapped_file_source> file_;
};

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
  std::string out;
};

void describe_bgmg_options(BgmgOptions& s) {
  LOG << "Options in effect (after applying default setting to non-specified parameters):";
  if (!s.bim.empty()) LOG << "\t--bim " << s.bim << " \\";
  if (!s.out.empty()) LOG << "\t--out " << s.out << " \\";
}

class BimFile {
public:
  BimFile(std::string filename) {
    LOG << "Reading " << filename << "...";

    std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
    std::stringstream ss;
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
      ss.str(str);
      ss >> chr_label >> snp >> gp >> bp >> a1 >> a2; 
      if (ss.fail()) {
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

    LOG << "Found " << chr_label_.size() << " variants.\n";
  }

private:
  std::vector<int> chr_label_;
  std::vector<std::string> snp_;
  std::vector<float> gp_;
  std::vector<int> bp_;
  std::vector<std::string> a1_;
  std::vector<std::string> a2_;
  friend class BimFileParser;
};

void fix_and_validate(BgmgOptions& bgmg_options, po::variables_map& vm) {
  // Validate --bim option
  if (bgmg_options.bim.empty()) throw std::invalid_argument(std::string("ERROR: --bim option must be specified"));
  if (!boost::filesystem::exists(bgmg_options.bim)) {
    std::stringstream ss; ss << "ERROR: input file " << bgmg_options.bim << " does not exist";
    throw std::runtime_error(ss.str());
  }

  // Validate --out option
  if (bgmg_options.out.empty())
    throw std::invalid_argument(std::string("ERROR: --out option must be specified"));
}

int main(int argc, char *argv[]) {
  try {
    BgmgOptions bgmg_options;
    po::options_description po_options("BGMG " VERSION " - Univariate and Bivariate causal mixture models for GWAS");
    po_options.add_options()
      ("help,h", "produce this help message")
      ("bim", po::value(&bgmg_options.bim), "Full path to .bim file that defines the reference set of SNPs")
      ("out", po::value(&bgmg_options.out)->default_value("bgmg"),
        "prefix of the output files; "
        "See README.md file for detailed description of file formats.")
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

      BimFile bim_file(bgmg_options.bim);

      auto analysis_finished = boost::posix_time::second_clock::local_time();
      LOG << "Analysis finished: " << analysis_finished;
      LOG << "Elapsed time: " << analysis_finished - analysis_started;
    }
    catch (std::exception& e) {
      LOG << e.what();
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
