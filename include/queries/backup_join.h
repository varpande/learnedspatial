#include "../Common.hpp"
#include "../boost/geometry.hpp"
// ---------------------------------------------------------------------------------------------------------------------
namespace bg = boost::geometry;
// ---------------------------------------------------------------------------------------------------------------------
//int crossing_number()
// ---------------------------------------------------------------------------------------------------------------------
/*
void join_query(const vector<Geography::Point> &points, const vector<Geography::Polygon> &polygons, int &count, int index) {
//	geos::geom::Point* p = pts.factory->createPoint(geos::geom::Coordinate(0.0, 0.0));
	for(uint64_t i = 0; i < points.size(); i++) {
//		geos::geom::Geometry* p = pts.factory->createPoint(geos::geom::Coordinate(points[i].y, points[i].x));
//		bg::model::point<double, 2, bg::cs::geographic<bg::degree>> p(points[i].y, points[i].x);
//		if(boost::geometry::within(p, polygon.polygon))
//		p->setX(points[i].y);
//		p->setY(points[i].x);
//		if(polygons.prep_geom.prepared_geometry[index]->contains(p))
//			count++;
	}
}
*/

/*
void join_query(const vector<Geography::Point> &points, const vector<Geography::Polygon> &polygons, int &count, int index) {
	std::cout << "Points size: " << points.size() << std::endl;
	for(uint64_t i = 0; i < points.size(); i++) {
		edge::range_type queryValue = points[i].y;
		itree<edge>::query_iterator iter = polygons[index].tree.qbegin(queryValue);
		itree<edge>::query_iterator end = polygons[index].tree.qend(queryValue);
		if(iter==end) continue;
		int cn = 0;
		while(iter != end) {
			int edge_id = iter->id();
			Geography::Point V = polygons[index].edges[edge_id];
			Geography::Point P = points[i];
			if ( ((polygons[index].edges[edge_id].y <= P.y) && (polygons[index].edges[edge_id+1].y > P.y)) || 
					 ((polygons[index].edges[edge_id].y > P.y) && (polygons[index].edges[edge_id+1].y <=  P.y))) {
				double vt = (double) (P.y  - polygons[index].edges[edge_id].y) / (polygons[index].edges[edge_id+1].y - polygons[index].edges[edge_id].y);

				if (P.x <  polygons[index].edges[edge_id].x + vt * (polygons[index].edges[edge_id+1].x - polygons[index].edges[edge_id].x))
					++cn;
			}
			iter++;
		}
		if (cn&1) count += 1;
	}
}

void join_query(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
//	std::cout << "Points size: " << points.size() << std::endl;
//	std::cout << "Point sample: [" << points[0].x << ", " << points[0].y << "]\n";
//	std::cout << "Geography::Polygon sample: [" << polygon.edges[0].x << ", " << polygon.edges[0].y << "]\n";
	for(uint64_t i = 0; i < points.size(); i++) {
		edge::range_type queryValue = points[i].y;
		itree<edge>::query_iterator iter = polygon.tree.qbegin(queryValue);
		itree<edge>::query_iterator end = polygon.tree.qend(queryValue);
		if(iter==end) continue;
		int cn = 0;
		while(iter != end) {
//			std::cout << "iter->low: " << iter->low() << ", iter->high: " << iter->high() << "\n";
			int edge_id = iter->id();
//			if(edge_id > polygon.edges.size() || edge_id < 0) break;
			Geography::Point V1 = polygon.edges[edge_id];
			Geography::Point V2 = polygon.edges[edge_id+1];
			Geography::Point P = points[i];

			if (((V1.y <= P.y) && (V2.y > P.y)) ||
         ((V1.y > P.y) && (V2.y <=  P.y))) {

				double vt = (double)(P.y  - V1.y) / (V2.y - V1.y);
				if (P.x <  V1.x + vt * (V2.x - V1.x))
					++cn;
			}
			iter++;
		}
		if (cn&1) count += 1;
	}
}

void full_join_query(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
  for(uint64_t i = 0; i < points.size(); i++) {
    edge::range_type queryValue = points[i].y;
//    itree<edge>::query_iterator iter = polygon.tree.begin();
//    itree<edge>::query_iterator end = polygon.tree.end();
		itree<edge>::const_iterator iter = polygon.tree.begin();
		itree<edge>::const_iterator end = polygon.tree.end();
    if(iter==end) continue;
    int cn = 0;
    while(iter != end) {
//      std::cout << "iter->low: " << iter->low() << ", iter->high: " << iter->high() << "\n";
      int edge_id = iter->id();
      if(edge_id > polygon.edges.size() || edge_id < 0) break;
      Geography::Point V1 = polygon.edges[edge_id];
      Geography::Point V2 = polygon.edges[edge_id+1];
      Geography::Point P = points[i];

      if (((V1.y <= P.y) && (V2.y > P.y)) ||
         ((V1.y > P.y) && (V2.y <=  P.y))) {

        double vt = (double)(P.y  - V1.y) / (V2.y - V1.y);
        if (P.x <  V1.x + vt * (V2.x - V1.x))
          ++cn;
      }
      iter++;
    }
    if (cn&1) count += 1;
  }
}
*/

void full_join_query_2(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
  for(uint64_t i = 0; i < points.size(); i++) {
		int cn = 0;
		Geography::Point P, V1, V2;
		P = points[i];
		V1 = polygon.edges[0];
		for(uint64_t j = 1; j < polygon.edges.size(); j++) {
			V2 = polygon.edges[j];

      if (((V1.y <= P.y) && (V2.y > P.y)) ||
         ((V1.y > P.y) && (V2.y <=  P.y))) {

        double vt = (double)(P.y  - V1.y) / (V2.y - V1.y);
        if (P.x <  V1.x + vt * (V2.x - V1.x))
          ++cn;
      }
			V1 = V2;
    }
    if (cn&1) count++;
  }
}

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
void join_umich(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
	for(uint64_t i = 0; i < points.size(); i++) {
		int counter = 0;
		Geography::Point p, p1, p2;
		p = points[i];
		p1 = polygon.edges[0];
		double xinters;
		for(uint64_t i = 1; i < polygon.edges.size(); i++) {
			p2 = polygon.edges[i];
			if (p.y > MIN(p1.y,p2.y)) {
				if (p.y <= MAX(p1.y,p2.y)) {
					if (p.x <= MAX(p1.x,p2.x)) {
						if (p1.y != p2.y) {
							xinters = (p.y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y)+p1.x;
							if (p1.x == p2.x || p.x <= xinters)
								counter++;
						}
					}
				}
			}
			p1 = p2;
		}
		if (counter % 2 != 0) count++;
	}
}

static bool inline intersect(double x1, double y1, double x2, double y2,
  						 				double x3, double y3, double x4, double y4) {
	if (x1 == x2) {
		return !(x3 == x4 && x1 != x3);
	} else if (x3 == x4) {
		return true;
	} else {
		double m1 = (y1-y2)/(x1-x2);
		double m2 = (y3-y4)/(x3-x4);
		return m1 != m2;
	}
}

static inline bool onSegment(Geography::Point &p, Geography::Point &q, Geography::Point &r) {
	if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
			q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
		return true;
	return false;
}

static inline int orientation(Geography::Point &p, Geography::Point &q, Geography::Point &r) {
	int val = (q.y - p.y) * (r.x - q.x) - 
						(q.x - p.x) * (r.y - q.y);
	if (val == 0) return 0;
	return (val > 0)? 1: 2;
}

static inline bool doIntersect(Geography::Point &p1, Geography::Point &q1, Geography::Point &p2, Geography::Point &q2) {
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1); 
  int o4 = orientation(p2, q2, q1);

	// General case 
  if (o1 != o2 && o3 != o4) 
        return true;
	
	// Special Cases 
  // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
  if (o1 == 0 && onSegment(p1, p2, q1)) return true; 
  
  // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
  if (o2 == 0 && onSegment(p1, q2, q1)) return true; 
  
  // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
  if (o3 == 0 && onSegment(p2, p1, q2)) return true; 
  
  // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
  if (o4 == 0 && onSegment(p2, q1, q2)) return true; 
  
  return false; // Doesn't fall in any of the above cases
}

void join_query_interval_tree(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
	// loop over points
  for(uint64_t i = 0; i < points.size(); i++) {

		int counter = 0;
		Geography::Point p, p1, p2, a,b,c,d;
		p = points[i];
		double xinters, x1, y1, x2, y2;

		// query the interval tree
		double query_value = p.y;
//		itree<edge>::query_iterator iter;
		itree<edge>::query_iterator iter = polygon.tree.qbegin(query_value);
		itree<edge>::query_iterator end = polygon.tree.qend(query_value);
		if(i < 10)
			std::cout << "query_value: " << query_value << ", intervals: ";

		while(iter != end) {
			x1 = iter->get_x1();
		  x2 = iter->get_x2();
			y1 = iter->get_y1();
      y2 = iter->get_y2();
			if(i < 10)
				std::cout << "[" << iter->low() << ", " << iter->high() << "], ";

/*
			if(p.y > MIN(y1 , y2)) {
				if (p.y <= MAX(y1, y2)) {
					if (p.x <= MAX(x1,x2)) {
						if (y1 != y2) {
							xinters = (p.y-y1)*(x2-x1)/(y2-y1)+x1;
							if (x1 == x2 || p.x <= xinters)
								counter++;
						}
					}
				}
			}
*/
			if((p.y > MIN(y1 , y2)) && (p.y <= MAX(y1, y2)) && (p.x <= MAX(x1,x2)) && (y1 != y2)) {
				xinters = (p.y-y1)*(x2-x1)/(y2-y1)+x1;
				if (x1 == x2 || p.x <= xinters)
					counter++;
			}


//			if(intersect(x1,y1,x2,y2,p.x,p.y,polygon.box.to.x,p.y)) counter++;
//			if(!intersect(x1,y1,x2,y2,p.x,p.y,90.0,p.y)) counter++;
			iter++;
/*
			p1 = iter->get_p1();
			p2 = iter->get_p2();
			if (p.y > MIN(p1.y,p2.y)) {
				if (p.y <= MAX(p1.y,p2.y)) {
					if (p.x <= MAX(p1.x,p2.x)) {
						if (p1.y != p2.y) {
							xinters = (p.y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y)+p1.x;
							if (p1.x == p2.x || p.x <= xinters)
								counter++;	
						}
					}
				}
			}
			iter++;
*/
		}
		if (counter%2) count++;
		if(i < 10)
			std::cout << std::endl;
	}

}

void test_sorted(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
	std::vector<pair<uint64_t, Geography::Point>> edges_copy;
	for(uint64_t i = 0; i < polygon.edges.size(); i++) {
		edges_copy.push_back(make_pair(i,polygon.edges[i]));
	}
	for(uint64_t i = 0; i < points.size(); i++) {
		int counter = 0;
		Geography::Point p = points[i];
	}
}

void new_interval(const vector<Geography::Point> &points, const Geography::Polygon &polygon, int &count) {
	Geography::Point p, p1, p2;
	vector<Interval<double, edge>> results;
	results.reserve(1e6);
	double xinters;
	for(uint64_t i = 0; i < points.size(); i++) {
		int counter = 0;
		p = points[i];
		results = polygon.tree.findOverlapping(p.y, p.y);
		for(uint64_t j = 0; j < results.size(); j++) {
			edge e = results[j].value;
			p1 = e.from;
      p2 = e.to;

			if((p.y > MIN(p1.y , p2.y)) && (p.y <= MAX(p1.y , p2.y)) &&
				 (p.x <= MAX(p1.x , p2.x)) && (p1.y != p2.y)) {
				xinters = (p.y - p1.y) * (p2.x - p1.x)/(p2.y - p1.y) + p1.x;
				if (p1.x == p2.x || p.x <= xinters)
					counter++;
			}
		}
		if(counter%2) count++;
	}
}
