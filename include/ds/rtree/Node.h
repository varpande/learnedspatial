#ifndef RTREE_INCLUDE_NODE_H_
#define RTREE_INCLUDE_NODE_H_

#include "Rectangle.h"

#include <vector>
#include <limits>
// ---------------------------------------------------------------------------------------------------------------------
namespace STRTree {
// ---------------------------------------------------------------------------------------------------------------------
template<uint32_t MAX_ENTRY_COUNT>
class Node {
 public:
  bool if_leaf;
  uint32_t entry_count;
  uint32_t bounds_count;
  std::array<Node<MAX_ENTRY_COUNT> *, MAX_ENTRY_COUNT> nodes; // Contains ids when is_leaf
  std::array<double, MAX_ENTRY_COUNT * 4> bounds;

  Node(int max_entry, bool leaf) {
    entry_count = 0;
    bounds_count = 0;
    if_leaf = leaf;
  }

  void inner_add_entry(Node *n, double minx, double miny, double maxx, double maxy) {
    nodes[entry_count] = n;

    bounds[bounds_count] = minx;
    bounds[bounds_count + 1] = miny;
    bounds[bounds_count + 2] = maxx;
    bounds[bounds_count + 3] = maxy;

    entry_count++;
    bounds_count += 4;
  }

  void leaf_add_entry(int id, double minx, double miny, double maxx, double maxy) {
    nodes[entry_count] = (Node<MAX_ENTRY_COUNT> *) (uint64_t) id;

    bounds[bounds_count] = minx;
    bounds[bounds_count + 1] = miny;
    bounds[bounds_count + 2] = maxx;
    bounds[bounds_count + 3] = maxy;

    entry_count++;
    bounds_count += 4;
  }

  bool is_leaf() { return (if_leaf); }

  double get_minx() {
    double minx = std::numeric_limits<double>::max();
    for (uint64_t i = 0; i < bounds_count; i = i + 4) {
      minx = std::min(minx, bounds[i]);
    }
    return minx;
  }

  double get_miny() {
    double miny = std::numeric_limits<double>::max();
    for (uint64_t i = 1; i < bounds_count; i = i + 4) {
      miny = std::min(miny, bounds[i]);
    }
    return miny;
  }

  double get_maxx() {
    double maxx = std::numeric_limits<double>::lowest();
    for (uint64_t i = 2; i < bounds_count; i = i + 4) {
      maxx = std::max(maxx, bounds[i]);
    }
    return maxx;
  }

  double get_maxy() {
    double maxy = std::numeric_limits<double>::lowest();
    for (uint64_t i = 3; i < bounds_count; i = i + 4) {
      maxy = std::max(maxy, bounds[i]);
    }
    return maxy;
  }

  void print_node() {
    std::cout << std::setprecision(15);
    std::cout << "===================================== Printing Node ======================================\n";
    std::cout << "entryCount: " << entry_count << "\n";
    double mbr_minx = get_minx();
    double mbr_miny = get_miny();
    double mbr_maxx = get_maxx();
    double mbr_maxy = get_maxy();
    std::cout << "Node bounds are: (" << mbr_minx << ", " << mbr_miny << ", " << mbr_maxx << ", " << mbr_maxy
              << ")\n";
    for (uint64_t i = 0; i < bounds_count; i = i + 4) {
      std::cout << "entry_num: " << i << " => (" << bounds[i] << ", " << bounds[i + 1] << ", " << bounds[i + 2]
                << ", " << bounds[i + 3] << ")\n";
    }
    std::cout << "=================================== Printing Node End ====================================\n";
  }
};
// ---------------------------------------------------------------------------------------------------------------------
} // namespace STRTree

#endif  // RTREE_INCLUDE_NODE_H_
