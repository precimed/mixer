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
#include <string>

#include <boost/program_options.hpp>
#include <boost/noncopyable.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/filesystem.hpp>
#include <boost/utility.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "bgmg.h"   // VERSION is defined here

namespace po = boost::program_options;

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

  BgmgCpp(int context_id) : context_id_(context_id) {}

  void init(std::string bim_file, std::string frq_file, std::string chr_labels, std::string trait1_file, std::string trait2_file, std::string exclude, std::string extract) {
    handle_errror(bgmg_init(context_id_, bim_file.c_str(), frq_file.c_str(), chr_labels.c_str(), trait1_file.c_str(), trait2_file.c_str(), exclude.c_str(), extract.c_str()));
  }

  void convert_plink_ld(std::string plink_ld_gz, std::string plink_ld_bin) {
    handle_errror(bgmg_convert_plink_ld(context_id_, plink_ld_gz.c_str(), plink_ld_bin.c_str()));
  }

private:
  void handle_errror(int error_code) {
    if (error_code < 0) throw std::runtime_error(bgmg_get_last_error());
  }

  int context_id_;
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
  std::string bfile;
  std::string bim;
  std::string frq;
  std::string chr_labels;
  std::string out;
  std::string plink_ld;
  std::string trait1;
  std::string exclude;
  std::string extract;
  float r2min;
  float ldscore_r2min;
};

void describe_bgmg_options(BgmgOptions& s) {
  LOG << "Options in effect (after applying default setting to non-specified parameters):";
  if (!s.bim.empty()) LOG << "\t--bim " << s.bim << " \\";
  if (!s.frq.empty()) LOG << "\t--frq " << s.frq << " \\";
  if (!s.out.empty()) LOG << "\t--out " << s.out << " \\";
  if (!s.plink_ld.empty()) LOG << "\t--plink-ld " << s.plink_ld << " \\";
  if (!s.trait1.empty()) LOG << "\t--trait1 " << s.trait1 << " \\";
  if (!s.exclude.empty()) LOG << "\t--exclude " << s.exclude << " \\";
  if (!s.extract.empty()) LOG << "\t--extract " << s.extract << " \\";
}

void fix_and_validate(BgmgOptions& bgmg_options, po::variables_map& vm) {
  // Validate --bim / --bim-chr option
  if (bgmg_options.bim.empty() && bgmg_options.bfile.empty())
    throw std::invalid_argument(std::string("ERROR: --bim or --bfile must be specified"));

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
  if (bgmg_options.frq.empty())
    throw std::invalid_argument(std::string("ERROR: --frq must be specified"));

  if (!bgmg_options.bfile.empty()) // ignore the remaining validation - this indicates that we built LD structure
    return;

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
      ("bfile", po::value(&bgmg_options.bfile), "Path to .fam/.bim/.bed file that defines the reference genotypes for LD structure estimatiion.")
      ("bim", po::value(&bgmg_options.bim), "Path to .bim file that defines the reference set of SNPs. Optionally, if input files are split per chromosome, use @ to specify location of chromosome label.")
      ("frq", po::value(&bgmg_options.frq), "Path to .frq file that defines the minor allele frequency for the reference set of SNPs. Optionally, if input files are split per chromosome, use @ to specify location of chromosome label.")
      ("plink-ld", po::value(&bgmg_options.plink_ld), "Path to plink .ld.gz file to convert into BGMG binary format.")
      ("chr-labels", po::value(&bgmg_options.chr_labels), "Set of chromosome labels. Defaults to '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22'")
      ("out", po::value(&bgmg_options.out)->default_value("bgmg"),
        "prefix of the output files; "
        "See README.md file for detailed description of file formats.")
        ("trait1", po::value(&bgmg_options.trait1), "Path to .sumstats.gz file for the trait to analyze")
      ("exclude", po::value(&bgmg_options.exclude)->default_value(""), "File with a set of SNP rs# to exclude from the analysis")
      ("extract", po::value(&bgmg_options.extract)->default_value(""), "File with a set of SNP rs# to use in the analysis; this is optional, by default use all available markers")
      ("r2min", po::value(&bgmg_options.r2min)->default_value(0.05), "Threshold for LD r2 estimation.")
      ("ldscore-r2min", po::value(&bgmg_options.ldscore_r2min)->default_value(0.05), "Threshold for LD scores in LD r2 estimation.")
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

      if (!bgmg_options.bfile.empty()) {
        bgmg_calc_ld_matrix(bgmg_options.bfile.c_str(),
                            bgmg_options.out.c_str(), bgmg_options.r2min, bgmg_options.ldscore_r2min, 0, 0);
      } else {
        const int context_id = 0;
        BgmgCpp bgmg_cpp_interface(context_id);
        bgmg_cpp_interface.init(bgmg_options.bim, bgmg_options.frq, bgmg_options.chr_labels, bgmg_options.trait1, std::string(), bgmg_options.exclude, bgmg_options.extract);

        if (!bgmg_options.plink_ld.empty()) {
          bgmg_cpp_interface.convert_plink_ld(bgmg_options.plink_ld, bgmg_options.out);
        }
      }

      auto analysis_finished = boost::posix_time::second_clock::local_time();
      LOG << "Analysis finished: " << analysis_finished;
      LOG << "Elapsed time: " << analysis_finished - analysis_started;
    }
    catch (std::exception& e) {
      std::cerr << "Exception  : " << e.what() << "\n";
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
