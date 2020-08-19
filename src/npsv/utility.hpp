#pragma once

#include <stdexcept>
#include <string>

inline void assert_throw(const bool cond, const std::string& text,
                         const std::string& file, const int line) {
  if (!cond) {
    throw std::runtime_error(text + ". In file: " + file +
                             " on line: " + std::to_string(line));
  }
}

#define pyassert(cond, text) assert_throw(cond, text, __FILE__, __LINE__)