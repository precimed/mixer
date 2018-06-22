#pragma once

#include "boost/utility.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

#include <fstream>
#include <vector>

#define LOG Logger::singleton()
#define MAX_LOG_LINES 10000000

class LoggerImpl : boost::noncopyable {
public:
  static LoggerImpl& singleton() {
    static LoggerImpl instance("bgmg.log");
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

private:
  int log_count_;  // ideally this should be boost::atomic<int>, but it's not a big deal if we loose some counts due to threading issues

  Logger(): log_count_() {
    (*this) << "============= new session =============";

  }
};
