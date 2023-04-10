#pragma once
// ---------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <cfloat>
#include <fcntl.h>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <chrono>
#include <thread>
#include <cassert>
#include <functional>
#include <algorithm>
#include <atomic>
#include <limits>

#include "boost/geometry.hpp"
// ---------------------------------------------------------------------------------------------------------------------
using namespace std;
// ---------------------------------------------------------------------------------------------------------------------
uint64_t GetFileLength(const string &file_name) {
  int fileFD = open(file_name.c_str(), O_RDONLY);
  if (fileFD < 0) {
    cerr << "Unable to open file: '" << file_name << "'" << endl; // You can check errno to see what happend
    throw;
  }
  if (fcntl(fileFD, F_GETFL) == -1) {
    cerr << "Unable to call fcntl on file: '" << file_name << "'" << endl; // You can check errno to see what happend
    throw;
  }
  struct stat st;
  fstat(fileFD, &st);
  close(fileFD);
  return st.st_size;
}

bool DoesFileExist(const std::string &file_name) {
  std::ifstream f(file_name.c_str());
  return f.good();
}

template<class POINT>
vector<POINT> LoadPoints(const string &file_name) {
  uint64_t length = GetFileLength(file_name);
  vector<POINT> points(length / sizeof(POINT));
  ifstream in(file_name);
  in.read((char *) &points[0], length);
  if (!in.good()) {
    cerr << "Could not read input file: '" << file_name << "'" << endl;
    throw;
  }
  return points;
}

template<class POINT>
void SavePoints(const string &file_name, const vector<POINT> &points) {
  uint64_t length = sizeof(POINT) * points.size();
  ofstream out(file_name);
  out.write((char *) &points[0], length);
  if (!out.good()) {
    cerr << "Could not write file: '" << file_name << "'" << endl;
    throw;
  }
}

inline vector<Geography::Rectangle> LoadRectanglesCsv(const string &file_name) {
#ifdef CACHE_CSV
  // Try using caching
 ifstream in_cache(string(CACHE_CSV) + file_name);
 if (in_cache.good()) {
    return LoadPoints<Geography::Rectangle>(string(CACHE_CSV) + file_name);
 }
#endif

  vector<Geography::Rectangle> result;
  result.reserve(1e6);
  ifstream in(file_name);
  if (!in.good()) {
    cerr << "Could not read input file: '" << file_name << "'" << endl;
    throw;
  }

  while (in.good()) {
    Geography::Rectangle rectangle;
    char comma;
    in >> rectangle.from.x >> comma >> rectangle.from.y >> comma >> rectangle.to.x >> comma >> rectangle.to.y;
    if (!in.good()) {
      break;
    }
    result.push_back(rectangle);
  }

#ifdef CACHE_CSV
  // Create cache
 SavePoints(string(CACHE_CSV) + file_name, result);
#endif

  return result;
}

inline vector<Geography::Circle> LoadDistanceCsv(const string &file_name) {
#ifdef CACHE_CSV
  // Try using caching
//   ifstream in_cache(string(CACHE_CSV) + file_name);
//   if (in_cache.good()) {
//      return LoadPoints<Geography::Rectangle>(string(CACHE_CSV) + file_name);
//   }
#endif

  vector<Geography::Circle> result;
  result.reserve(1e6);
  ifstream in(file_name);
  if (!in.good()) {
    cerr << "Could not read input file: '" << file_name << "'" << endl;
    throw;
  }

  while (in.good()) {
    Geography::Circle circle;
    char comma;
    in >> circle.center.x >> comma >> circle.center.y >> comma >> circle.radius;
    if (!in.good()) {
      break;
    }
    result.push_back(circle);
  }

#ifdef CACHE_CSV
  // Create cache
//   SavePoints(string(CACHE_CSV) + file_name, result);
#endif

  return result;
}

static void parse_csv(std::ifstream &csv_input,
                      std::function<void(std::vector<std::string> &)> callback) {

  using namespace boost;
  typedef tokenizer<escaped_list_separator<char>> Tokenizer;

  std::vector<std::string> tokens;
  std::string line;

  while (getline(csv_input, line)) {
    Tokenizer tok(line);
    tokens.assign(tok.begin(), tok.end());
    callback(tokens);
  }
}

static void read_wkt(const std::string &wkt, Geography::Polygon &p, Geography::Rectangle &r) {
  const char *begin = wkt.c_str() + 9;
  char *next = nullptr;
  double longitude, latitude;

  double minx = DBL_MAX;
  double miny = DBL_MAX;
  double maxx = -DBL_MAX;
  double maxy = -DBL_MAX;

  for (auto iter = begin; iter != begin + wkt.length(); iter++) {
    longitude = strtod(iter, &next);
    iter = ++next;
    latitude = strtod(iter, &next);

    if (latitude < minx) minx = latitude;
    if (latitude > maxx) maxx = latitude;
    if (longitude < miny) miny = longitude;
    if (longitude > maxy) maxy = longitude;

    Geography::Point single_point;
    single_point.x = latitude;
    single_point.y = longitude;
    p.edges.push_back(single_point);

    if (*next == ')') break;
    else iter = next;
  }

  r.from.x = minx;
  r.from.y = miny;
  r.to.x = maxx;
  r.to.y = maxy;
}

static void parse_wkt_polygons(vector<Geography::Polygon> &polygons,
                               const std::string &file) {
  std::ifstream in(file);
  if (!in.is_open()) {
    std::cerr << "unable to open " << file << std::endl;
    exit(1);
  }

  bool osm = (file.find("countries") != std::string::npos);
  if (osm) std::cout << "OSM is set" << std::endl;
  parse_csv(in, [&](std::vector<std::string> &fields) {
    Geography::Polygon p;
    Geography::Rectangle r;

    if (osm)
      read_wkt(fields[0], p, r);
    else
      read_wkt(fields[2], p, r);
    p.box = r;
    std::cout.precision(cout.precision(std::numeric_limits<double>::max_digits10));
//    std::cout << polygons.size() << ", Geography::Rectangle[" << r.from.x << ", " << r.from.y << ", " << r.to.x << ", " << r.to.y << "]" << std::endl;

    polygons.push_back(p);
  });
}
// ---------------------------------------------------------------------------------------------------------------------
