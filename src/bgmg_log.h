#pragma once

#include "boost/utility.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/algorithm/string.hpp"
#include <fstream>
#include <vector>
#include <chrono>

#define LOG Logger::singleton()
#define MAX_LOG_LINES 10000000
#define BGMG_THROW_EXCEPTION(x) throw (x)

// A timer that fire an event each X milliseconds.
class SimpleTimer {
public:
  SimpleTimer(int period_ms) : start_(std::chrono::system_clock::now()), period_ms_(period_ms) {}

  int elapsed_ms() {
    auto delta = (std::chrono::system_clock::now() - start_);
    auto delta_ms = std::chrono::duration_cast<std::chrono::milliseconds>(delta);
    return delta_ms.count();
  }

  bool fire() {
    if (elapsed_ms() < period_ms_) return false;
    start_ = std::chrono::system_clock::now();
    return true;
  }
private:
  std::chrono::time_point<std::chrono::system_clock> start_;
  int period_ms_;
};

class LoggerImpl : boost::noncopyable {
public:
  static LoggerImpl& singleton() {
    static LoggerImpl instance;
    return instance;
  }

  static LoggerImpl& singleton_void() {
    static LoggerImpl instance;
    return instance;
  }

private:
  LoggerImpl() : log_file_() {}
  LoggerImpl(std::string path) : log_file_() {
    try {
      init(path);
    } catch (...) {
      // ignore errors (maybe two instances are trying to write to the same bgmg.log file?)
    }

    if (!log_file_.is_open()) {
      std::cerr << "unable to initialize logging to " << path;
    }
  }

public:
  bool is_initialized() {
    return log_file_.is_open();
  }

  void init(std::string path) {
    if (log_file_.is_open()) {
      log_file_.close();
    }

    log_file_.open(path, std::ios::app);

    if (log_file_.is_open()) {
      boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%Y%m%d %H:%M:%S.%f");
      log_file_.imbue(std::locale(log_file_.getloc(), facet));
    }
    // boost::posix_time::time_facet *facet = new boost::posix_time::time_facet("%d-%b-%Y %H:%M:%S");
    // std::cout.imbue(std::locale(std::cout.getloc(), facet));
    // log_file_.imbue(std::locale(log_file_.getloc(), facet));
    // std::cout << std::setprecision(9);
  }

  template <typename T>
  LoggerImpl& operator<< (const T& rhs) {
    // std::cout << rhs;
    if (log_file_.is_open()) {
      log_file_ << rhs;
      log_file_.flush();
    }
    return *this;
  }

  template <typename T>
  LoggerImpl& operator<< (const std::vector<T>& rhs) {
    for (int i = 0; i < rhs.size(); i++)
      (*this) << ((i>0) ? " " : "") << rhs[i];
    return (*this);
  }

private:
  std::ofstream log_file_;
};

class Logger : boost::noncopyable {
public:
  static Logger& singleton() {
    static Logger instance;
    return instance;
  }

  template <typename T>
  LoggerImpl& operator<< (const T& rhs) {
    auto now = boost::posix_time::microsec_clock::local_time();

    int log_count = log_count_++;
    if (log_count == MAX_LOG_LINES) {
      LoggerImpl::singleton() << "\n" << now << "\t" << "Too many lines written to the log; log throtling enabled.\n";
    }

    if (log_count >= MAX_LOG_LINES) {
      return LoggerImpl::singleton_void();
    }

    return (LoggerImpl::singleton() << "\n" << now << "\t" << rhs);
  }

  static std::vector<std::string> tokenize_message(std::string message) {
    std::vector<std::string> tokens;
    const std::string separators = "\n\r";
    boost::trim_if(message, boost::is_any_of(separators));
    boost::split(tokens, message, boost::is_any_of(separators), boost::token_compress_on);
    return tokens;
  }

private:
  int log_count_;  // ideally this should be boost::atomic<int>, but it's not a big deal if we loose some counts due to threading issues

  Logger(): log_count_() {
    (*this) << "============= new session =============";

  }
};

