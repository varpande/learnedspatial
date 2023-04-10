#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"

#include "queries/distance.h"
#include "queries/join.h"
#include "partitioning_techniques/FixedGrid.hpp"
#include "partitioning_techniques/AdaptiveGrid.hpp"
#include "partitioning_techniques/KdTreePartitioning.hpp"
#include "partitioning_techniques/QuadtreePartitioning.hpp"
#include "partitioning_techniques/STRPartitioning.hpp"
#include "partition_cells/BinarySearchY.hpp"
#include "partition_cells/Spline.hpp"

#include <unordered_map>
#include <map>
#include <iostream>
#include <algorithm>
// ---------------------------------------------------------------------------------------------------------------------
#define VALIDATE_MATERIALIZE 1
#define VALIDATE_COUNT 2
#define BENCHMARK_MATERIALIZE 3
#define BENCHMARK_COUNT 4
#define BENCHMARK_DISTANCE 5
#define BENCHMARK_POINT 6
#define BENCHMARK_JOIN 7
#ifndef BENCHMARK_TYPE
#define BENCHMARK_TYPE BENCHMARK_COUNT
#endif
#if BENCHMARK_TYPE == VALIDATE_MATERIALIZE
#define RUN_FUNC ValidateMaterialize
#endif
#if BENCHMARK_TYPE == VALIDATE_COUNT
#define RUN_FUNC ValidateCount
#endif
#if BENCHMARK_TYPE == BENCHMARK_MATERIALIZE
#define RUN_FUNC BenchmarkMaterialize
#endif
#if BENCHMARK_TYPE == BENCHMARK_COUNT
#define RUN_FUNC BenchmarkCount
#endif
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
#define RUN_FUNC BenchmarkDistance
#endif
#if BENCHMARK_TYPE == BENCHMARK_POINT
#define RUN_FUNC BenchmarkPoint
#endif
#if BENCHMARK_TYPE == BENCHMARK_JOIN
#define RUN_FUNC BenchmarkJoin
#endif
#ifndef RUN_FUNC
#error Invalid BENCHMARK_TYPE, please use: { "VALIDATE_MATERIALIZE", "VALIDATE_COUNT", "BENCHMARK_MATERIALIZE", "BENCHMARK_COUNT", "BENCHMARK_DISTANCE", "BENCHMARK_POINT", "BENCHMARK_JOIN"}
#endif
// ---------------------------------------------------------------------------------------------------------------------
#ifndef PARTITION_SIZE
#define PARTITION_SIZE 1000
#endif
#ifndef PARTITION_SIZE_GRIDS
#define PARTITION_SIZE_GRIDS 1000
#endif
// ---------------------------------------------------------------------------------------------------------------------
string g_points_file;
string g_rectangles_file;
// ---------------------------------------------------------------------------------------------------------------------
// Fallback for weird bugs ..
struct FullScanEngine {
  const vector<Geography::Point> &points;

  FullScanEngine(const vector<Geography::Point> &points, const Geography::Rectangle &)
      : points(points) {
  }

  void LookUp(const Geography::Rectangle &rectangle, vector<Geography::Point> &result) const {
    assert(result.empty());
    for (const Geography::Point &p : points) {
      if (rectangle.from.x <= p.x && p.x <= rectangle.to.x && rectangle.from.y <= p.y && p.y <= rectangle.to.y) {
        result.push_back(p);
      }
    }
  }

  uint32_t Count(const Geography::Rectangle &rectangle) const {
    uint32_t result = 0;
    for (const Geography::Point &p : points) {
      if (rectangle.from.x <= p.x && p.x <= rectangle.to.x && rectangle.from.y <= p.y && p.y <= rectangle.to.y) {
        result++;
      }
    }
    return result;
  }
};
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void BenchmarkMaterialize(const vector<Geography::Point> &points,
                          const vector<Geography::Rectangle> &rectangles,
                          const Geography::Rectangle &bounding_box) {
  COMPETITOR competitor;
  uint64_t build_ns = RunWithTime([&]() {
    competitor.Build(points, bounding_box);
  });

  vector<Geography::Point> result;
  result.reserve(points.size() * 0.01);
  uint64_t count = 0;

  uint64_t run_ns = RunWithTime([&]() {
    for (uint64_t i = 0; i < rectangles.size(); i++) {
      competitor.LookUp(rectangles[i], result);
      count += result.size();
      result.clear();
    }
  });
  //@formatter:off
   cout << "RES " << competitor.GetName()
        << " COUNT: " << count
        << " points_file: " << g_points_file
        << " rectangles_file: " << g_rectangles_file
        << " #points: " << points.size()
        << " #lookups: " << rectangles.size()
        << " memory(byte): " << competitor.GetUsedMemory()
        << " build(ms): " << (build_ns / 1000 / 1000)
        << " ns/lookup: " << run_ns / rectangles.size() << endl;
   //@formatter:on
}
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void BenchmarkCount(const vector<Geography::Point> &points,
                    const vector<Geography::Rectangle> &rectangles,
                    const Geography::Rectangle &bounding_box) {
  COMPETITOR competitor;

#ifdef PRINT_STATS
    number_of_cells = 0;
    wrongly_scanned_points_count = 0;
    actual_scanned_points_count = 0;
    cells_intersected = 0;
#endif

  uint64_t build_ns = RunWithTime([&]() {
    competitor.Build(points, bounding_box);
  });

  uint64_t count = 0;

  uint64_t run_ns = RunWithTime([&]() {
    for (uint64_t i = 0; i < rectangles.size(); i++) {
      count += competitor.Count(rectangles[i]);
    }
  });

//	double total_scanned_points_count = ((double)(actual_scanned_points_count + wrongly_scanned_points_count) / (double) rectangles.size());
  double result_count = (((double) count) / (double) rectangles.size());
//	double scan_overhead = total_scanned_points_count / result_count;
//	double scan_overhead = ((double)(actual_scanned_points_count + wrongly_scanned_points_count) / (double) rectangles.size());
//	double scan_overhead = ((double)(actual_scanned_points_count / (double) rectangles.size()) / (double)result_count);
//	double avg_points = ((double) points.size()) / (double) number_of_cells;

//	double cells = ((double) cells_intersected) / (double) rectangles.size();

#ifdef PRINT_STATS
  double total_scanned_points_count = ((double)(actual_scanned_points_count + wrongly_scanned_points_count) / (double) rectangles.size());
  double scan_overhead = ((double)(actual_scanned_points_count / (double) rectangles.size()) / (double)result_count);
  double cells = ((double) cells_intersected) / (double) rectangles.size();
#endif

  //@formatter:off
   cout << "RES " << competitor.GetName()
        << " COUNT: " << count
        << " points_file: " << g_points_file
        << " rectangles_file: " << g_rectangles_file
        << " #points: " << points.size()
        << " #lookups: " << rectangles.size()
        << " memory(byte): " << competitor.GetUsedMemory()
        << " build(ms): " << (build_ns / 1000 / 1000)
        << std::fixed
        << " result_count: " << result_count
//		  << " average_points: " << avg_points
        << " ns/lookup: " << run_ns / rectangles.size()
#if PRINT_STATS
        << " num_cells: " << number_of_cells
        << " actual_scanned_points_count: " << actual_scanned_points_count
        << " wrongly_scanned_points_count: " << wrongly_scanned_points_count
        << " total_scanned_points_count: " << total_scanned_points_count
        << " scan_overhead: " << scan_overhead
        << " cells_intersected: " << cells
#endif
		  << endl;

//	cout << "wrongly_scanned_points_count:   " << wrongly_scanned_points_count / rectangles.size() << endl;
//	cout << "actual_scanned_points_count:    " << actual_scanned_points_count / rectangles.size() << endl;
//	cout << "total_scanned_points_count:     " << (actual_scanned_points_count + wrongly_scanned_points_count) / rectangles.size() << endl;

	cout << "------------------------------------------------------------------------\n";
   //@formatter:on
}
// ---------------------------------------------------------------------------------------------------------------------
void CheckResult(vector<Geography::Point> &result, vector<Geography::Point> &reference_result) {
  if (result.size() != reference_result.size()) {
    cout << "ERROR: different result sizes:" << endl;
    cout << "Expected: " << reference_result.size() << endl;
    cout << "Actual: " << result.size() << endl;
    throw;
  }
  sort(result.begin(), result.end());
  sort(reference_result.begin(), reference_result.end());
  for (uint32_t i = 0; i < result.size(); i++) {
    if (result[i] != reference_result[i]) {
      cout << "ERROR: different result content at [" << i << "]" << endl;
      throw;
    }
  }
}
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void ValidateMaterialize(const vector<Geography::Point> &points,
                         const vector<Geography::Rectangle> &rectangles,
                         const Geography::Rectangle &bounding_box) {
  COMPETITOR competitor;
  competitor.Build(points, bounding_box);

  // Compare everything against binary search
  {
    FixedGridPartitioning<BinarySearchYCell, 1000> reference;
    reference.Build(points, bounding_box);
    for (uint64_t i = 0; i < rectangles.size(); i++) {
      vector<Geography::Point> result;
      vector<Geography::Point> reference_result;
      competitor.LookUp(rectangles[i], result);
      reference.LookUp(rectangles[i], reference_result);
      CheckResult(result, reference_result);
    }
  }

  cout << COMPETITOR::GetName() << ": materialize test ok!" << endl;
}
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void ValidateCount(const vector<Geography::Point> &points,
                   const vector<Geography::Rectangle> &rectangles,
                   const Geography::Rectangle &bounding_box) {
  COMPETITOR competitor;
  competitor.Build(points, bounding_box);

  // Compare everything against binary search
  {
    FixedGridPartitioning<BinarySearchYCell, 1000> reference;
    reference.Build(points, bounding_box);
    for (uint64_t i = 0; i < rectangles.size(); i++) {
      uint32_t result = competitor.Count(rectangles[i]);
      uint32_t reference_result = reference.Count(rectangles[i]);
      if (result != reference_result) {
        cout << "Got different results: reference: " << reference_result << " result: " << result << endl;
        throw;
      }
    }
  }

  cout << COMPETITOR::GetName() << ": count test ok!" << endl;
}
// ---------------------------------------------------------------------------------------------------------------------
void GetAllTheStatistics(const vector<Geography::Point> &points) {
  unordered_map<double, uint64_t> x_duplicate_map;
  unordered_map<double, uint64_t> y_duplicate_map;
  map<Geography::Point, uint64_t> points_duplicate_map;
  uint64_t x_duplicate_count = 0;
  uint64_t y_duplicate_count = 0;
  uint64_t point_duplicate_count = 0;

  for (auto &iter : points) {
    if (x_duplicate_map.count(iter.x)) {
      x_duplicate_count++;
    }
    x_duplicate_map[iter.x]++;

    if (y_duplicate_map.count(iter.y)) {
      y_duplicate_count++;
    }
    y_duplicate_map[iter.y]++;

    if (points_duplicate_map.count(iter)) {
      point_duplicate_count++;
    }
    points_duplicate_map[iter]++;
  }

  cout << "point_duplicate_count: " << point_duplicate_count << endl;
  cout << "x_duplicate_count: " << x_duplicate_count << endl;
  cout << "y_duplicate_count: " << y_duplicate_count << endl;
  cout << "unique_location: " << points_duplicate_map.size() << endl;
  cout << "x_unique_location: " << x_duplicate_map.size() << endl;
  cout << "y_unique_location: " << y_duplicate_map.size() << endl;

  // vector<pair<Geography::Point, uint64_t>> histogram;
  // for(auto iter : points_duplicate_map) {
  //    histogram.push_back(iter);
  // }
  // sort(histogram.begin(), histogram.end(), [](const pair<Geography::Point, uint64_t>& lhs, const pair<Geography::Point, uint64_t>& rhs){return lhs.second>rhs.second;});
}
// ---------------------------------------------------------------------------------------------------------------------
void FilterPoints(vector<Geography::Point> &points) {
#ifdef POINTS_FILTER
  assert(POINTS_FILTER>0);
  uint32_t write_idx = 0;
  for (uint32_t read_idx = 0; read_idx<points.size(); read_idx++) {
     if (read_idx % POINTS_FILTER == 0) {
        points[write_idx] = points[read_idx];
        write_idx++;
     }
  }
  points.resize(points.size() / POINTS_FILTER);
#endif
  (void) points;
}
// ---------------------------------------------------------------------------------------------------------------------
void FilterQueryPoints(vector<Geography::Point> &points) {
#ifdef QUERY_POINTS_FILTER
  assert(QUERY_POINTS_FILTER>0);
  uint32_t write_idx = 0;
  for (uint32_t read_idx = 0; read_idx<points.size(); read_idx++) {
     if (read_idx % QUERY_POINTS_FILTER == 0) {
        points[write_idx] = points[read_idx];
        write_idx++;
     }
  }
  points.resize(points.size() / QUERY_POINTS_FILTER);
#endif
  (void) points;
}

// ---------------------------------------------------------------------------------------------------------------------
void FilterRectangles(vector<Geography::Rectangle> &rectangles) {
#ifdef RECTANGLES_FILTER
  assert(RECTANGLES_FILTER>0);
  uint32_t write_idx = 0;
  for (uint32_t read_idx = 0; read_idx<rectangles.size(); read_idx++) {
     if (read_idx % RECTANGLES_FILTER == 0) {
        rectangles[write_idx] = rectangles[read_idx];
        write_idx++;
     }
  }
  rectangles.resize(rectangles.size() / RECTANGLES_FILTER);
#endif
  (void) rectangles;
}
// ---------------------------------------------------------------------------------------------------------------------
void FilterCircles(vector<Geography::Circle> &query_points, vector<Geography::Rectangle> &rectangles) {
#ifdef CIRCLE_FILTER
  assert(CIRCLE_FILTER>0);
  uint32_t write_idx = 0;
  for (uint32_t read_idx = 0; read_idx<rectangles.size(); read_idx++) {
     if (read_idx % CIRCLE_FILTER == 0) {
        rectangles[write_idx] = rectangles[read_idx];
                query_points[write_idx] = query_points[read_idx];
        write_idx++;
     }
  }
  rectangles.resize(rectangles.size() / CIRCLE_FILTER);
    query_points.resize(query_points.size() / CIRCLE_FILTER);
#endif
  (void) rectangles;
  (void) query_points;
}
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void BenchmarkDistance(const vector<Geography::Point> &points,
                       const vector<Geography::Circle> &query_points,
                       const vector<Geography::Rectangle> &query_rectangles,
                       const Geography::Rectangle &bounding_box) {
  COMPETITOR competitor;
  uint64_t build_ns = RunWithTime([&]() {
    competitor.Build(points, bounding_box);
  });

  vector<Geography::Point> rectangle_result;
  rectangle_result.reserve(points.size() * 0.01);
  vector<Geography::Point> distance_result;
  distance_result.reserve(points.size() * 0.01);
  uint64_t count = 0;
  uint64_t time_refine = 0;
  uint64_t time_filter = 0;

  uint64_t run_ns = RunWithTime([&]() {
    for (uint64_t i = 0; i < query_rectangles.size(); i++) {
      auto start = chrono::high_resolution_clock::now();
      // filter
      competitor.LookUp(query_rectangles[i], rectangle_result);
      auto finish = chrono::high_resolution_clock::now();
      time_filter += chrono::duration_cast<chrono::nanoseconds>(finish - start).count();

      // refine
      start = chrono::high_resolution_clock::now();
      distance_query(rectangle_result, query_points[i].center, query_points[i].radius, distance_result);
      finish = chrono::high_resolution_clock::now();
      time_refine += chrono::duration_cast<chrono::nanoseconds>(finish - start).count();

      // Bookkeeping
      count += distance_result.size();
      distance_result.clear();
      rectangle_result.clear();
    }
  });
  //@formatter:off
    cout << "RES " << competitor.GetName()
         << " COUNT: " << count
         << " points_file: " << g_points_file
         << " distances_file: " << g_rectangles_file
         << " #points: " << points.size()
         << " #lookups: " << query_rectangles.size()
         << " #avg_result: " << (count / query_rectangles.size())
         << " memory(byte): " << competitor.GetUsedMemory()
         << " build(ms): " << (build_ns / 1000 / 1000)
         << " total_time(ns): " << run_ns
         << " filter_time(ns): " << time_filter
         << " refine_time(ns): " << time_refine
         << " ns/lookup: " << run_ns / query_rectangles.size()
         << " ns/refine: " << time_refine / query_rectangles.size() << endl;
    //@formatter:on
}
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void BenchmarkPoint(const vector<Geography::Point> &points,
                    const vector<Geography::Point> &rectangles,
                    const Geography::Rectangle &bounding_box) {
  COMPETITOR competitor;
  uint64_t build_ns = RunWithTime([&]() {
    competitor.Build(points, bounding_box);
  });

  vector<Geography::Point> result;
  result.reserve(points.size() * 0.01);
  uint64_t count = 0;

  uint64_t run_ns = RunWithTime([&]() {
    for (uint64_t i = 0; i < rectangles.size(); i++) {
      Geography::Rectangle r;
      r.from = rectangles[i];
      r.to = rectangles[i];
      competitor.PointLookUp(rectangles[i], result);
      count += result.size();
      result.clear();
    }
  });
  //@formatter:off
   cout << "RES " << competitor.GetName()
        << " COUNT: " << count
        << " points_file: " << g_points_file
        << " point_queries_file: " << g_rectangles_file
        << " #points: " << points.size()
        << " #lookups: " << rectangles.size()
        << " memory(byte): " << competitor.GetUsedMemory()
        << " build(ms): " << (build_ns / 1000 / 1000)
        << " ns/lookup: " << run_ns / rectangles.size() << endl;
   //@formatter:on
}
// ---------------------------------------------------------------------------------------------------------------------
template<class COMPETITOR>
void BenchmarkJoin(const vector<Geography::Point> &points,
                   vector<Geography::Polygon> &polygons,
                   const Geography::Rectangle &bounding_box) {
  typedef std::numeric_limits<double> dbl;
  COMPETITOR competitor;
  uint64_t build_ns = RunWithTime([&]() {
    competitor.Build(points, bounding_box);
  });

  vector<Geography::Point> rectangle_result;
  rectangle_result.reserve(points.size() * 0.01);
  vector<int> join_result(polygons.size(), 0);
  uint64_t count = 0;
  uint64_t time_refine = 0;

  uint64_t run_ns_interval = RunWithTime([&]() {
    // iterate over the polygons, use the bounding box to filter
    // and the refine using interval tree of each polygon
    for (uint64_t i = 0; i < polygons.size(); i++) {
      // filter
      competitor.LookUp(polygons[i].box, rectangle_result);

      // refine
      auto start = chrono::high_resolution_clock::now();
      points_in_polygon_join(rectangle_result, polygons[i], join_result[i]);
      auto finish = chrono::high_resolution_clock::now();
      time_refine += chrono::duration_cast<chrono::nanoseconds>(finish - start).count();

      // bookkeeping
      count += join_result[i];
      rectangle_result.clear();
    }
  });

  //@formatter:off
    cout << "RES " << competitor.GetName()
         << " COUNT: " << count
         << " points_file: " << g_points_file
         << " polygon_file: " << g_rectangles_file
         << " #points: " << points.size()
         << " #polygons: " << polygons.size()
         << " #avg_result: " << (count / polygons.size())
         << " memory(byte): " << competitor.GetUsedMemory()
         << " build(ms): " << (build_ns / 1000 / 1000)
         << " total_time(ns): " << run_ns_interval
         << " refine_time(ns): " << time_refine
         << " total_time(s): " << run_ns_interval / (1e9)
         << " refine_time(s): " << time_refine / 1e9
         << " ns/lookup: " << run_ns_interval / polygons.size() << endl;
    //@formatter:on

  return;
}
// ---------------------------------------------------------------------------------------------------------------------
auto PrintConfig(int argc, char **argv) {
  cout << "\nLoading:" << endl;
  cout << "========" << endl;

  string type;
  ostringstream warnings;
#if BENCHMARK_TYPE == VALIDATE_MATERIALIZE
  type =  "Validate Materialize :)";
#ifdef NDEBUG
 warnings << "WARNING: NDEBUG is defined in validation run!!!" << endl;
#endif
#endif
#if BENCHMARK_TYPE == VALIDATE_COUNT
  type = "Validate Count :)";
#ifdef NDEBUG
 warnings << "WARNING: NDEBUG is defined in validation run!!!" << endl;
#endif
#endif
#if BENCHMARK_TYPE == BENCHMARK_MATERIALIZE
  type = "Benchmark Materialize o.O";
#ifndef NDEBUG
 warnings << "WARNING: NDEBUG is not defined in benchmark run!!!" << endl;
#endif
#endif
#if BENCHMARK_TYPE == BENCHMARK_COUNT
  type = "Benchmark Count o.O";
#ifndef NDEBUG
  warnings << "WARNING: NDEBUG is not defined in benchmark run!!!" << endl;
#endif
#endif
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  type = "Benchmark Distance o.O";
#ifndef NDEBUG
  warnings << "WARNING: NDEBUG is not defined in benchmark run!!!" << endl;
#endif
#endif
#if BENCHMARK_TYPE == BENCHMARK_POINT
  type = "Benchmark Point Query o.O";
#ifndef NDEBUG
  warnings << "WARNING: NDEBUG is not defined in benchmark run!!!" << endl;
#endif
#endif
#if BENCHMARK_TYPE == BENCHMARK_Join
  type = "Benchmark Join o.O";
#ifndef NDEBUG
  warnings << "WARNING: NDEBUG is not defined in benchmark run!!!" << endl;
#endif
#endif
#ifdef POINTS_FILTER
  warnings << "WARNING: Using POINTS_FILTER '" << POINTS_FILTER << "'" << endl;
#endif
#ifdef RECTANGLES_FILTER
  warnings << "WARNING: Using RECTANGLES_FILTER '" << RECTANGLES_FILTER << "'" << endl;
#endif
#ifdef CACHE_CSV
  warnings << "WARNING: Using CACHE_CSV '" << CACHE_CSV << "'" << endl;
#endif
#ifdef PRINT_STATS
  warnings << "WARNING: Printing stats is enabled, the benchmark run time measurements may be noisy and will incur extra overhead." << endl;
#endif

  if (argc != 3) {
    cout << "usage: " << argv[0] << " <points_file> <rectangles_file>" << endl;
    exit(-1);
  }

  g_points_file = argv[1];
  g_rectangles_file = argv[2];

  // Points
  vector<Geography::Point> points;
  cout << left << setw(50) << setfill('.') << ("Points: ") << " " << flush;
  bool has_sorted_file = false;
  uint64_t points_load_time = RunWithTime([&]() {
    if (DoesFileExist(g_points_file + ".sorted_x")) {
      points = LoadPoints<Geography::Point>(g_points_file + ".sorted_x");
      has_sorted_file = true;
    } else {
      points = LoadPoints<Geography::Point>(g_points_file);
    }
  });


  // Sorting
  cout << left << setw(50) << setfill('.') << "Sorting: " << " " << flush;
  uint64_t sort_time = RunWithTime([&]() {
    if (!has_sorted_file) {
      SortPointsX(points);
      SavePoints(g_points_file + ".sorted_x", points);
    }
  });
  cout << (sort_time / 1000 / 1000 / 1000.0) << endl;

  // Filter
  uint64_t points_filter_time = 0;
#ifdef POINTS_FILTER
  cout << left << setw(50) << setfill('.') << ("Filter: ") << " " << flush;
 points_filter_time = RunWithTime([&]() {
    FilterPoints(points);
 });
 cout << (points_filter_time / 1000 / 1000 / 1000.0) << endl;
#endif

#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  vector<Geography::Point> points_r(points.size());
  cout << left << setw(50) << setfill('.') << ("Generated Points Radians: ") << " " << flush;
  uint64_t points_r_load_time = RunWithTime([&]() {
    if (DoesFileExist(g_points_file + ".radians")) {
      points_r = LoadPoints<Geography::Point>(g_points_file + ".radians");
    } else {
      points_r = ConvertToRadian(points);
      SavePoints(g_points_file + ".radians", points_r);
      points.clear();
    }
#ifdef POINTS_FILTER
    FilterPoints(points_r);
#endif
  });
  cout << (points_r_load_time / 1000 / 1000 / 1000.0) << endl;
#endif

#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  vector<Geography::Circle> query_points;
  vector<Geography::Rectangle> distance_rectangles;
  cout << left << setw(50) << setfill('.') << ("Distances: ") << " " << flush;
  uint64_t rectangles_load_time = RunWithTime([&]() {
    query_points = LoadDistanceCsv(g_rectangles_file);
  });
  cout << (rectangles_load_time / 1000 / 1000 / 1000.0) << endl;

  cout << left << setw(50) << setfill('.') << ("Distances Convert: ") << " " << flush;
  vector<Geography::Circle> query_points_r;
  uint64_t rectangles_convert_time = RunWithTime([&]() {
    GetRectanglesFromDistance(query_points, distance_rectangles);
          FilterCircles(query_points, distance_rectangles);
    query_points_r = ConvertToRadian(query_points);
    query_points.clear();
  });
  cout << (rectangles_convert_time / 1000 / 1000 / 1000.0) << endl;
#elif BENCHMARK_TYPE == BENCHMARK_POINT
  vector<Geography::Point> query_points;
  cout << left << setw(50) << setfill('.') << ("Query Points: ") << " " << flush;
  uint64_t rectangles_load_time = RunWithTime([&]() {
          query_points = LoadPoints<Geography::Point>(g_rectangles_file);
          FilterQueryPoints(query_points);
  });
  cout << (rectangles_load_time / 1000 / 1000 / 1000.0) << endl;
#elif BENCHMARK_TYPE == BENCHMARK_JOIN
  vector<Geography::Polygon> polygons;
  cout << left << setw(50) << setfill('.') << ("Polygons: ") << " " << flush;
  uint64_t rectangles_load_time = RunWithTime([&]() {
    parse_wkt_polygons(polygons, g_rectangles_file);
    for(uint64_t i = 0; i < polygons.size(); i++) {
      polygons[i].construct_tree();
    }
  });
  cout << (rectangles_load_time / 1000 / 1000 / 1000.0) << endl;
#else
  // Rectangles
  cout << left << setw(50) << setfill('.') << ("Rectangles: ") << " " << flush;
  vector<Geography::Rectangle> rectangles;
  uint64_t rectangles_load_time = RunWithTime([&]() {
    rectangles = LoadRectanglesCsv(g_rectangles_file);
    FilterRectangles(rectangles);
  });
  cout << (rectangles_load_time / 1000 / 1000 / 1000.0) << endl;
#endif
  // Bounding box
  cout << left << setw(50) << setfill('.') << "Bounding box: " << " " << flush;
  Geography::Rectangle bounding_box;
  uint64_t bounding_box_time = RunWithTime([&]() {
    bounding_box = GetBoundingBox(points);
  });
  cout << (bounding_box_time / 1000 / 1000 / 1000.0) << endl;

  // Total
  cout << left << setw(50) << setfill('.') << "Total setup time(s): " << " " << flush;
  cout << ((points_load_time + sort_time + points_filter_time + rectangles_load_time + bounding_box_time) / 1000 / 1000
      / 1000.0) << endl;

  cout << "\nConfig:" << endl;
  cout << "=======" << endl;
  cout << "Type:          " << type << endl;
  if (warnings.str().size() == 0) {
    cout << "Warnings:      none" << endl;
  } else {
    cout << "\033[1;31m" << warnings.str() << "\033[0m";
  }
  cout << "Points file:   '" << g_points_file << "'" << endl;
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  cout << "Point.size():  " << points_r.size() << endl;
#else
  cout << "Point.size():  " << points.size() << endl;
#endif
  cout << "Rectangles file:   '" << g_rectangles_file << "'" << endl;
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  cout << "Distance.size(): " << distance_rectangles.size() << endl;
#elif BENCHMARK_TYPE == BENCHMARK_POINT
  cout << "Query_Points.size(): " << query_points.size() << endl;
#elif BENCHMARK_TYPE == BENCHMARK_JOIN
  cout << "Polygons.size(): " << polygons.size() << endl;
#else
  cout << "Rectangles.size(): " << rectangles.size() << endl;
#endif
#ifdef POINTS_FILTER
  cout << "POINTS_FILTER: " << POINTS_FILTER << endl;
#endif
#ifdef RECTANGLES_FILTER
  cout << "RECTANGLES_FILTER: " << RECTANGLES_FILTER << endl;
#endif

  cout << "\nRun:" << endl;
  cout << "====" << endl;
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  return make_tuple(points_r, query_points_r, distance_rectangles, bounding_box);
#elif BENCHMARK_TYPE == BENCHMARK_POINT
  return make_tuple(points, query_points, bounding_box);
#elif BENCHMARK_TYPE == BENCHMARK_JOIN
  return make_tuple(points, polygons, bounding_box);
#else
  return make_tuple(points, rectangles, bounding_box);
#endif
}
// ---------------------------------------------------------------------------------------------------------------------
int main(int args, char **argv) {

  // Query Specific Input
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  auto[points, query_points, query_rectangles, bounding_box] = PrintConfig(args, argv);
#elif BENCHMARK_TYPE == BENCHMARK_POINT
  auto[points, query_points, bounding_box] = PrintConfig(args, argv);
#elif BENCHMARK_TYPE == BENCHMARK_JOIN
  auto[points, polygons, bounding_box] = PrintConfig(args, argv);
#else
  auto [points, rectangles, bounding_box] = PrintConfig(args, argv);
#endif

  // Run Benchmarks
  // Only FixedGrid is activated by default
#if BENCHMARK_TYPE == BENCHMARK_DISTANCE
  RUN_FUNC<FixedGridPartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<KdTreePartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<STRPartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<QuadtreePartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);

  RUN_FUNC<FixedGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<KdTreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<STRPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);
  RUN_FUNC<QuadtreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, query_rectangles, bounding_box);

#elif BENCHMARK_TYPE == BENCHMARK_POINT
  RUN_FUNC<FixedGridPartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<KdTreePartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<STRPartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<QuadtreePartitioning<SplineCell, PARTITION_SIZE>>(points, query_points, bounding_box);

  RUN_FUNC<FixedGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<KdTreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<QuadtreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, bounding_box);
  RUN_FUNC<STRPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, query_points, bounding_box);

#elif BENCHMARK_TYPE == BENCHMARK_JOIN
  RUN_FUNC<FixedGridPartitioning<SplineCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<SplineCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<KdTreePartitioning<SplineCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<STRPartitioning<SplineCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<QuadtreePartitioning<SplineCell, PARTITION_SIZE>>(points, polygons, bounding_box);

  RUN_FUNC<FixedGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<KdTreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<STRPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, polygons, bounding_box);
  RUN_FUNC<QuadtreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, polygons, bounding_box);

// else run range query benchmarks
#else
  RUN_FUNC<FixedGridPartitioning<SplineCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<SplineCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<KdTreePartitioning<SplineCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<STRPartitioning<SplineCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<QuadtreePartitioning<SplineCell, PARTITION_SIZE>>(points, rectangles, bounding_box);

  RUN_FUNC<FixedGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<AdaptiveGridPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<KdTreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<STRPartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
  RUN_FUNC<QuadtreePartitioning<BinarySearchYCell, PARTITION_SIZE>>(points, rectangles, bounding_box);
#endif

  return 0;
}
// ---------------------------------------------------------------------------------------------------------------------
