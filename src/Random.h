/*
 * Random number generator based on Mersenne Twister
 * Random.h
 * Released under the MIT License, MMC051 Contributor
 */
#pragma once
#include <random>

namespace Random {
  inline std::mt19937& engine() {
      static thread_local std::mt19937 eng([]() {
          std::random_device rd;
          return std::mt19937(rd());
      }());
      return eng;
  }

  inline int get(int min, int max) {
      std::uniform_int_distribution<int> dist(min, max);
      return dist(engine());
  }
}
