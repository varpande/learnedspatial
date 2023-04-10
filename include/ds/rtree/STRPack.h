#ifndef RTREE_INCLUDE_STRPACK_H_
#define RTREE_INCLUDE_STRPACK_H_

#include "ds/rtree/Data.h"
#include "ds/rtree/Node.h"
#include "ds/rtree/Rectangle.h"

#include <cmath>
#include <vector>
#include <queue>
#include <algorithm>

//#include "schema.h"
// ---------------------------------------------------------------------------------------------------------------------
// this class implements the STR packing algorithm for R-tree as presented in the following
// publication:
// "STR: A simple and efficient algorithm for R-tree packing."
// Leutenegger, Scott T., Mario A. Lopez, and Jeffrey Edgington.
// Proceedings 13th international conference on data engineering. IEEE, 1997.
// Algorithmic credits belong to the original authors.

namespace STRTree {
// ---------------------------------------------------------------------------------------------------------------------
template<uint32_t MAX_ENTRY_COUNT_IN_NODE>
class STRPack {
 public:
  typedef STRTree::Node<MAX_ENTRY_COUNT_IN_NODE> Node;

  bool built;
  int max_node_entries;
  int min_node_entries;

  std::vector<Data> data;
  Node *root_node;

  STRPack(int max_entries, int min_entries) {
    if (max_entries != MAX_ENTRY_COUNT_IN_NODE) {
      throw;
    }

    if (max_entries < 2) {
      max_entries = 50;
    }

    if (min_entries < 1 || min_entries > max_entries / 2) {
      min_entries = max_entries / 2;
    }

    init(max_entries, min_entries);
  }

  STRPack() {
    init(10, 5);
  }

  void init(int max_entries, int min_entries) {
    built = false;
    max_node_entries = max_entries;
    min_node_entries = min_entries;
    root_node = nullptr;
  }

  static bool compare_on_x(Data a, Data b) {
    double a_center_x = a.object.center_x();
    double b_center_x = b.object.center_x();
    return (a_center_x < b_center_x);
  }

  static bool compare_on_y(Data a, Data b) {
    double a_center_y = a.object.center_y();
    double b_center_y = b.object.center_y();
    return (a_center_y < b_center_y);
  }

  void add(const Rectangle &r, int id) {
    Data new_entry;
    new_entry.object = r;
    new_entry.object_id = id;
    data.push_back(new_entry);
  }

  void build() {
    if (!built) {
      build_original_str();
      built = true;
    }
  }

  // this function builds the tree following guidelines from the publication
  void build_original_str() {
    // sort all entries on the x-dimension
    std::sort(data.begin(), data.end(), compare_on_x);

    // p is number of leaf level pages according to the original publication
    int p = (int) std::ceil(data.size() / ((double) max_node_entries));

    // s (or num_slices) is the number of vertical slices
    int num_slices = (int) std::ceil(std::sqrt(p));

    int partition_size = num_slices * max_node_entries;

    int levels = compute_levels() + 1;

    int num_partitions = data.size() / partition_size;

    // sort each partition on the y-dimension
    for (int i = 0; i < (num_partitions - 1) * partition_size; i = i + partition_size) {
      auto begin = data.begin() + i;
      auto end = begin + partition_size;
      std::sort(begin, end, compare_on_y);
    }
    // sort the last partition on the y-dimension
    std::sort(data.begin() + ((num_partitions - 1) * partition_size), data.end(), compare_on_y);

    build_tree(levels);
    built = true;
  }

  // this is a variant that only sort on the x-dimension
  void build_variant_str() {
    std::sort(data.begin(), data.end(), compare_on_x);
    int levels = compute_levels() + 1;
    build_tree(levels);
  }

  int compute_levels() {
    double base = max_node_entries;
    double x = data.size();
    return (int) (std::log(x) / std::log(base));
  }

  void build_tree(int num_levels) {
    std::vector<std::vector<Node *>> nodes(num_levels);
    build_leaves(nodes);

    // Leaves are level 1, after leaves have been built we build the inner nodes
    // starting from level 2
    for (int i = 2; i <= num_levels; i++) {
      int loop_size = nodes[i - 2].size();
      Node *new_node = nullptr;
      Node *n = nullptr;

      // loop size determines how many nodes are at a level below the current level
      // we iterate over those nodes below forming new nodes per max_node_entries
      // if loop_size is <= max_node_entries then we know this will be the last level => root
      for (int j = 0; j < loop_size; j++) {

        // previous node is now full, make a new node
        if ((j % max_node_entries) == 0) {
          new_node = new Node(max_node_entries, false);
          nodes[i - 1].push_back(new_node);
          if (loop_size <= max_node_entries)
            root_node = new_node;
        }

        n = nodes[i - 2][j];
        double minx = n->get_minx();
        double miny = n->get_miny();
        double maxx = n->get_maxx();
        double maxy = n->get_maxy();
        new_node->inner_add_entry(n, minx, miny, maxx, maxy);
      }
    }
  }

  void build_leaves(std::vector<std::vector<Node *>> &nodes) {
    int loop_size = std::floor(data.size() / max_node_entries);
    int remaining_elements = data.size() % max_node_entries;

    // make loop_size number of nodes
    for (int i = 0; i < loop_size; i++) {
      Node *new_node = new Node(max_node_entries, true);
      for (int j = 0; j < max_node_entries; j++) {
        int index = (i * max_node_entries) + j;
        Rectangle r = data[index].object;
        int index_id = data[index].object_id;
        new_node->leaf_add_entry(index_id, r.minx, r.miny, r.maxx, r.maxy);
      }
      nodes[0].push_back(new_node);
    }

    // make last node with remaining_elements number of elements
    int begin = (loop_size * max_node_entries);
    int end = data.size();

    if (begin != end) {
      Node *last_leaf_node = new Node(max_node_entries, true);
      for (int i = begin; i < end; i++) {
        Rectangle r = data[i].object;
        int index_id = data[i].object_id;
        last_leaf_node->leaf_add_entry(index_id, r.minx, r.miny, r.maxx, r.maxy);
      }
      nodes[0].push_back(last_leaf_node);
    }
  }

  // query functions follow
  void query_intersects(const Rectangle &r, std::vector<int> &result) {
    intersects(r, result, root_node);
  }

  void intersects(const Rectangle &r, std::vector<int> &result, Node *n) {
    for (int i = 0; i < n->bounds_count; i = i + 4) {
      if (Rectangle::intersects(r.minx, r.miny, r.maxx, r.maxy, n->bounds[i], n->bounds[i + 1], n->bounds[i + 2],
                                n->bounds[i + 3])) {
        if (n->is_leaf()) {
          result.push_back((uint64_t) n->nodes[i >> 2]);
        } else {
          intersects(r, result, n->nodes[i >> 2]);
        }
      }
    }
  }

  void intersects_queue(Rectangle r, std::vector<int> &result) {
    std::queue<Node *> q;
    Node *rootNode = root_node;
    q.push(rootNode);
    while (q.empty() == false) {
      Node *n = q.front();
      for (int i = 0; i < n->bounds_count; i = i + 4) {
        if (Rectangle::intersects(r.minx, r.miny, r.maxx, r.maxy, n->bounds[i], n->bounds[i + 1],
                                  n->bounds[i + 2], n->bounds[i + 3])) {
          if (n->is_leaf()) {
            result.push_back((uint64_t) n->nodes[i >> 2]);
          } else
            q.push(n->nodes[i >> 2]);
        }
      }
      q.pop();
    }
  }
};
// ---------------------------------------------------------------------------------------------------------------------
} // namespace STRTree

#endif  // RTREE_INCLUDE_STRPACK_H_
