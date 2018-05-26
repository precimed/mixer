#pragma once

#include "boost/utility.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"

#include <fstream>
#include <vector>

#define LOG Logger::singleton()

class LoggerImpl : boost::noncopyable {
public:
  static LoggerImpl& singleton() {
    static LoggerImpl instance("bgmg.log");
    return instance;
  }
private:
  LoggerImpl(std::string path) : log_file_(path, std::ios::app) {
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%Y%m%d %H:%M:%S.%f");
    log_file_.imbue(std::locale(log_file_.getloc(), facet));
    // boost::posix_time::time_facet *facet = new boost::posix_time::time_facet("%d-%b-%Y %H:%M:%S");
    // std::cout.imbue(std::locale(std::cout.getloc(), facet));
    // log_file_.imbue(std::locale(log_file_.getloc(), facet));
    // std::cout << std::setprecision(9);
  }
public:
  template <typename T>
  LoggerImpl& operator<< (const T& rhs) {
    // std::cout << rhs;
    log_file_ << rhs;
    log_file_.flush();
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
    return (LoggerImpl::singleton() << "\n" << now << "\t" << rhs);
  }

private:
  Logger() {
    (*this) << "============= new session =============";

  }
};
