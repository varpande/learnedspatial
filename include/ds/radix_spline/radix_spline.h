#pragma once

#include "util.h"
#include "spline_util.h"

#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>

namespace spline {

template<class KeyType, class ValueType>
class RadixSpline {

 public:
  void Build(const std::vector<std::pair<KeyType, ValueType>> &data) {
    uint64_t fitting_coord_count = data.size() / 300; // XXX tuning parameter

    data_ = data;
    cdf = spline::buildCdf(data);

    // Compute the spline
    spline = spline::compressFunc(cdf, fitting_coord_count);
    spline_size = spline.size();

    // Store the spline x-coordinates in a radix index
    buildRadix();
  }

  // Wiki: https://en.wikipedia.org/wiki/Infimum_and_supremum
  // Returns the index to the greatest element that is less than or equal to the search key
  // (Excpetion: When the search key is smaller than the smallest element, 0 is returned)
  uint64_t GetIndexOfInfimum(const KeyType lookup_key) const {
    if (lookup_key <= data_[0].first) {
      return 0;
    }

    uint64_t estimate = segmentInterpolation(process(lookup_key), lookup_key);
    return util::linear_search(data_, lookup_key, estimate);
  }

  std::size_t size() const {
    return sizeof(*this) + spline.size() * sizeof(Coord) // for spline
        + ((1ull << num_radix_bits_) + 1) * sizeof(uint32_t) // for radix_hint
        + data_.size() * sizeof(KeyValue<KeyType>);
  }

  // Choose the most appropriate num_radix_bits_
  // Note: 20 for osm in SOSD
  void SetTuning(uint32_t num_radix_bits) {
    num_radix_bits_ = num_radix_bits;
  }

 private:
  std::vector<Coord> cdf;
  std::vector<Coord> spline;
  uint64_t spline_size;

  static constexpr uint32_t defaultValueRadixBits = 20;
  uint32_t num_radix_bits_ = 20;
  uint64_t n_;
  KeyType min_;
  KeyType max_;
  KeyType shift_bits_;
  uint32_t *radix_hint_;

  // Copy of data.
  std::vector<std::pair<KeyType, ValueType>> data_;

  inline uint64_t shift_bits(const uint64_t val) {
    const uint32_t clz = __builtin_clzl(val);
    if ((64 - clz) < num_radix_bits_)
      return 0;
    else
      return 64 - num_radix_bits_ - clz;
  }

  inline uint32_t shift_bits(const uint32_t val) {
    const uint32_t clz = __builtin_clz(val);
    if ((32 - clz) < num_radix_bits_)
      return 0;
    else
      return 32 - num_radix_bits_ - clz;
  }

  void buildRadix() {
    assert(num_radix_bits_);

    // Alloc the memory for the hints
    radix_hint_ = new uint32_t[(1ull << num_radix_bits_) + 1];

    // Compute the number of bits to shift with
    n_ = spline.size();
    min_ = spline.front().first;
    max_ = spline.back().first;
    shift_bits_ = shift_bits(max_ - min_);

    // Compute the hints
    radix_hint_[0] = 0;
    uint64_t prev_prefix = 0;
    for (uint64_t i = 0; i < n_; ++i) {
      uint64_t tmp = static_cast<uint64_t>(spline[i].first);
      uint64_t curr_prefix = (tmp - min_) >> shift_bits_;
      if (curr_prefix != prev_prefix) {
        for (uint64_t j = prev_prefix + 1; j <= curr_prefix; ++j)
          radix_hint_[j] = i;
        prev_prefix = curr_prefix;
      }
    }

    // Margin hint values
    for (; prev_prefix < (1ull << num_radix_bits_); ++prev_prefix)
      radix_hint_[prev_prefix + 1] = n_;
  }

  inline uint32_t process(uint64_t x) const {
    // Compute index.
    uint32_t index;
    const uint64_t p = (x - min_) >> shift_bits_;
    uint32_t begin = radix_hint_[p];
    uint32_t end = radix_hint_[p + 1];

    switch (end - begin) {
      case 0: index = end;
        break;
      case 1: index = (spline[begin].first >= x) ? begin : end;
        break;
      case 2: index = ((spline[begin].first >= x) ? begin : ((spline[begin + 1].first >= x) ? (begin + 1) : end));
        break;
      case 3: index =
                  ((spline[begin].first >= x) ? begin : ((spline[begin + 1].first >= x) ? (begin + 1) : ((spline[begin
                      + 2].first > x) ? (begin + 2) : end)));
        break;
      default:
        index = std::lower_bound(spline.begin() + begin,
                                 spline.begin() + end,
                                 x,
                                 [](const Coord &a, const uint64_t lookup_key) {
                                   return a.first < lookup_key;
                                 }) - spline.begin();
        break;
    }
    return index - 1;
  }

  double segmentInterpolation(uint64_t segment, const double x) const
  // get f(x) at segment
  {
    Coord down = spline[segment], up = spline[segment + 1];
    double slope = (down.second - up.second) / (down.first - up.first);
    return down.second + (x - down.first) * slope;
  }
};

}

