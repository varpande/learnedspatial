#include "radix_spline.h"
#include <map>
#include <functional>
#include <random>
// ---------------------------------------------------------------------------------------------------------------------
using namespace std;
// ---------------------------------------------------------------------------------------------------------------------
struct Point {
  double x, y;
};
// ---------------------------------------------------------------------------------------------------------------------
uint64_t LookUpInReference(map<uint64_t, Point> &reference, uint64_t key) {
  auto range = reference.equal_range(key);
  if (range.first == reference.end()) {
    return reference.size() - 1;
  }
  if (range.first == reference.begin()) {
    return 0;
  } else {
    if (range.first->first == key) {
      return distance(reference.begin(), range.first);
    } else {
      return distance(reference.begin(), range.first) - 1;
    }
  }
}
// ---------------------------------------------------------------------------------------------------------------------
uint64_t LookUpInReference(multimap<uint64_t, Point> &reference, uint64_t key) {
  auto range = reference.equal_range(key);
  if (range.first == reference.end()) {
    return reference.size() - 1;
  }
  if (range.first == reference.begin()) {
    return 0;
  } else {
    if (range.first->first == key) {
      return distance(reference.begin(), range.first);
    } else {
      return distance(reference.begin(), range.first) - 1;
    }
  }
}
// ---------------------------------------------------------------------------------------------------------------------
void WhiteBoardTets() {
  spline::RadixSpline<uint64_t, Point> radix_spline;

  vector<pair<uint64_t, Point>> data;
  data.push_back(make_pair(2, Point{0.0, 0.0}));
  data.push_back(make_pair(5, Point{0.0, 0.0}));
  data.push_back(make_pair(6, Point{0.0, 0.0}));
  data.push_back(make_pair(8, Point{0.0, 0.0}));
  radix_spline.Build(data);

  map<uint64_t, Point> reference;
  reference.insert(make_pair(2, Point{0.0, 0.0}));
  reference.insert(make_pair(5, Point{0.0, 0.0}));
  reference.insert(make_pair(6, Point{0.0, 0.0}));
  reference.insert(make_pair(8, Point{0.0, 0.0}));

  //   cout << "5 |-> " << LookUpInReference(reference, 5) << endl;
  //   cout << "4 |-> " << LookUpInReference(reference, 4) << endl;
  //   cout << "1 |-> " << LookUpInReference(reference, 1) << endl;

  //   cout << "5 |-> " << radix_spline.GetIndexOfInfimum(5) << endl;
  //   cout << "4 |-> " << radix_spline.GetIndexOfInfimum(4) << endl;
  //   cout << "1 |-> " << radix_spline.GetIndexOfInfimum(1) << endl;

  assert(LookUpInReference(reference, 5) == 1);
  assert(LookUpInReference(reference, 4) == 0);
  assert(LookUpInReference(reference, 1) == 0);

  assert(radix_spline.GetIndexOfInfimum(5) == 1);
  assert(radix_spline.GetIndexOfInfimum(4) == 0);
  assert(radix_spline.GetIndexOfInfimum(1) == 0);

  cout << "White board tests: passed :)" << endl;
}
// ---------------------------------------------------------------------------------------------------------------------
void RandomTest(int duplicate_level) {
  assert(duplicate_level == 0 || duplicate_level == 1 || duplicate_level == 2);

  multimap<uint64_t, Point> reference;

  vector<pair<uint64_t, Point>> data;
  spline::RadixSpline<uint64_t, Point> radix_spline;

  std::default_random_engine random_engine;
  random_engine.seed(8128);
  std::uniform_int_distribution<uint64_t> dist(0, 1000000);

  for (uint32_t i = 0; i < 10000; i++) {
    uint64_t key = dist(random_engine);
    if (duplicate_level == 2) {
      key /= 10000;
    }
    pair<uint64_t, Point> entry(key, Point{i * 1000.0, i * 1000.0});
    if (duplicate_level == 0 && reference.count(entry.first) > 0) {
      continue;
    }
    data.push_back(entry);
    reference.insert(entry);
  }
  sort(data.begin(),
       data.end(),
       [](const pair<uint64_t, Point> &lhs, const pair<uint64_t, Point> &rhs) { return lhs.first < rhs.first; });
  radix_spline.Build(data);

  for (uint32_t i = 0; i < 10000; i++) {
    uint64_t lookup_key = dist(random_engine);
    uint64_t reference_result = LookUpInReference(reference, lookup_key);
    uint64_t radix_result = radix_spline.GetIndexOfInfimum(lookup_key);

    if (reference_result != radix_result) {
      cout << "ref: " << reference_result << endl;
      cout << "radix: " << radix_result << endl;
    }
    assert(reference_result == radix_result);
  }

  cout << "Random tests (" << duplicate_level << "): passed :)" << endl;
}
// ---------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {
  WhiteBoardTets();
  RandomTest(0);
  RandomTest(1);
  RandomTest(2);

  return 0;
}
// ---------------------------------------------------------------------------------------------------------------------
