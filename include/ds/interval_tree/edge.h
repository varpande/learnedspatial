#pragma once
// -------------------------------------------------------------------------------------
#include "ds/geography/DataTypes.hpp"
#include "utils/IO.hpp"
#include "utils/Utils.hpp"
struct Geography::Point;
class edge {
 public:
  typedef double range_type;

  edge(range_type rHigh, range_type rLow, uint64_t sID, Geography::Point p1, Geography::Point p2)
      : m_rHigh(rHigh), m_rLow(rLow), m_sID(sID), p1(p1), p2(p2) {
    std::cout << "Inserting, low: " << rLow << ", high: " << rHigh << ", edge_id: " << sID << std::endl;
  }
  double low() const { return m_rLow; }
  double high() const { return m_rHigh; }
  uint64_t id() const { return m_sID; }
  Geography::Point get_p1() { return p1; }
  Geography::Point get_p2() { return p2; }

 private:
  range_type m_rHigh, m_rLow;
  uint64_t m_sID;
  Geography::Point p1;
  Geography::Point p2;
};
