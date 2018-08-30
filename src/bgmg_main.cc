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

// Perform concurrent parsing of a single text file, keeping track of lines.
// ResultT - result type
// ParserT - type that implements the actual parsing, exposed as follows:
//  void push(std::string& line);  - remember a line (don't parse it yet)
//  void parse();                  - do the actual parsing
//  void flush(ResultT* result, int first_line_index);   - append result to container
// ParserT must be reusable (e.i. support multiple sequences of push()->parse()->flush()->push()->...
// parse_input_in_parallel guaranties that each parse() will operate on a subsequent set of lines from the input file.
// 0-based index of the first line per chunk is passed to flush() as first_line_index.
template <typename ResultT, typename ParserT>
void parse_input_in_parallel(std::string filename, int lines_per_chunk, ResultT* result) {
  std::mutex read_access;
  std::mutex write_access;

  int global_line_no = 0;
  ifstream_or_cin stream_or_cin(filename);
  std::istream& infile = stream_or_cin.get_stream();

  int num_threads;
#pragma omp parallel
#pragma omp master  
  num_threads = omp_get_num_threads();

  auto func = [&infile, &read_access, &write_access, &global_line_no, lines_per_chunk, filename, result]() {
    int first_line_index = 0;
    ParserT parser;
    while (true) {
      {
        std::lock_guard<std::mutex> guard(read_access);
        if (infile.eof()) break;
        first_line_index = global_line_no;
        for (int line_index = 0; line_index < lines_per_chunk; line_index++) {
          std::string line;
          std::getline(infile, line);
          if (infile.eof()) break;
          global_line_no++;
          parser.push(line);
        }
      }

      parser.parse();

      {
        std::lock_guard<std::mutex> guard(write_access);
        parser.flush(result, first_line_index);
      }
    }
  };

  // The func may throw an exception if input file is malformed.
  // This exception will be re-thrown on the main thread.
  // http://stackoverflow.com/questions/14222899/exception-propagation-and-stdfuture
  std::vector<std::shared_future<void>> tasks;
  for (int i = 0; i < num_threads; i++) tasks.push_back(std::move(std::async(std::launch::async, func)));
  for (int i = 0; i < num_threads; i++) tasks[i].get();
}

class BimFile;
class BimFileParser {
public:
  void push(std::string& line) { lines_.push_back(line); }
  void parse() {
    const std::string separators = " \t\n\r";
    std::vector<std::string> tokens;

    chr_label_.resize(lines_.size());
    snp_.resize(lines_.size());
    gp_.resize(lines_.size());
    bp_.resize(lines_.size());
    a1_.resize(lines_.size());
    a2_.resize(lines_.size());

    for (int line_index = 0; line_index < lines_.size(); line_index++) {
      boost::trim_if(lines_[line_index], boost::is_any_of(separators));
      boost::split(tokens, lines_[line_index], boost::is_any_of(separators), boost::token_compress_on);

      try {
        chr_label_[line_index] = stoi(tokens[0]);  // must not use stream or boost::lexical_cast (result in some lock contention)
        snp_[line_index] = tokens[1];
        gp_[line_index] = stof(tokens[2]);
        bp_[line_index] = stoi(tokens[3]);
        a1_[line_index] = tokens[4];
        a2_[line_index] = tokens[5];
      }
      catch (...) {
        std::stringstream ss;
        ss << "Error parsing line #" << (line_index + 1) << ": " << lines_[line_index];
        if (line_index == 0) ss << "\nNB! bim files must not contain header line (that's the definition from plink).";
        throw std::runtime_error(ss.str());
      }
    }
  }
  void flush(BimFile* bimfile, int line_start_index);

private:
  std::string filename_;
  std::vector<std::string> lines_;
  std::vector<int> chr_label_;
  std::vector<std::string> snp_;
  std::vector<float> gp_;
  std::vector<int> bp_;
  std::vector<std::string> a1_;
  std::vector<std::string> a2_;
};

class BimFile {
public:
  BimFile(std::string filename) {
    LOG << "Reading " << filename << "...";
    parse_input_in_parallel<BimFile, BimFileParser>(filename, 100000, this);
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

void BimFileParser::flush(BimFile* bimfile, int first_line_index) {
  const int capacity = std::max(bimfile->chr_label_.size(), first_line_index + lines_.size());
  bimfile->chr_label_.resize(capacity);
  bimfile->snp_.resize(capacity);
  bimfile->gp_.resize(capacity);
  bimfile->bp_.resize(capacity);
  bimfile->a1_.resize(capacity);
  bimfile->a2_.resize(capacity);
  for (int i = 0; i < lines_.size(); i++) {
    bimfile->chr_label_[first_line_index + i] = chr_label_[i];
    bimfile->snp_[first_line_index + i] = snp_[i];
    bimfile->gp_[first_line_index + i] = gp_[i];
    bimfile->bp_[first_line_index + i] = bp_[i];
    bimfile->a1_[first_line_index + i] = a1_[i];
    bimfile->a2_[first_line_index + i] = a2_[i];
  }

  lines_.resize(0); // never reduce the actual capacity
}

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
