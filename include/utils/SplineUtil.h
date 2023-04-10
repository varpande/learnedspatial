#pragma once
// ---------------------------------------------------------------------------------------------------------------------
namespace spline {
// ---------------------------------------------------------------------------------------------------------------------
using segCoord = std::pair<double, uint32_t>;
using Coord = std::pair<double, double>;
using Errors = std::pair<double, double>; // could also be changed to support relative error

static double interpolate(const std::vector<Coord> &spline, double pos)
// evaluate the linear function found on spline for "pos"
{
  auto iter = lower_bound(spline.begin(), spline.end(), pos, [](const Coord &a, double b) { return a.first < b; });
  if (iter->first == pos)
    return iter->second;

  // Compute slope
  double dx = (iter + 0)->first - (iter - 1)->first;
  double dy = (iter + 0)->second - (iter - 1)->second;

  // Compute offset for the linear function
  double ofs = pos - (iter - 1)->first;
  return (iter - 1)->second + ofs * (dy / dx);
}

static Errors computeErrors(const std::vector<Coord> &cdf, const std::vector<Coord> &spline)
// returns the average and maximum errors
{
  // Encoding: errs.first -> avg, errs.second -> max
  Errors errs;
  for (auto elem : cdf) {
    double pos = elem.first;
    double estimate = interpolate(spline, pos);
    double real = elem.second;
    double error = estimate - real;
    if (error < 0)
      error = -error;
    if (error > errs.second)
      errs.second = error;
    errs.first += error;
  }
  errs.first /= cdf.size();
  return errs;
}

static void printErrors(const std::vector<Coord> &cdf, const std::vector<Coord> &spline) __attribute__((unused));
static void printErrors(const std::vector<Coord> &cdf, const std::vector<Coord> &spline)
// print the errors of the spline to the cdf
{
  Errors errs = computeErrors(cdf, spline);
  std::cout << "cdf_size: " << cdf.size() << std::endl << "spline_size: " << spline.size() << std::endl
            << "Errors for spline - (average, max): (" << errs.first << ", " << errs.second << ")" << std::endl;
}

// Compare derivations. -1 if (a->b) is more shallow than (a->c), +1 if it is more steeper
static int cmpDevs(const Coord &a, const Coord &b, const Coord &c) {
  double dx1 = b.first - a.first, dx2 = c.first - a.first;
  double dy1 = b.second - a.second, dy2 = c.second - a.second;
  if ((dy1 * dx2) < (dy2 * dx1)) // == (dy1/dx1) < (dy2/dx2)
    return -1;
  if ((dy1 * dx2) > (dy2 * dx1))
    return 1;
  return 0;
}

// lowerLimit and upperLimit span a corridor, if the iter point is not in it we add it
// Otherwise, it is kind of on the same line and we dont need to add it
// On each iteration the corridor is narrowed
// Epsilon specifies the error bound: larger epsilon -> larger error corridor
static std::vector<Coord> tautString(const std::vector<Coord> &data, double /*maxValue*/, double epsilon) {
  // Add the first point
  std::vector<Coord> result;
  if (data.empty())
    return result;
  result.push_back(data[0]);
  if (data.size() == 1)
    return result;

  // Taut the string
  auto iter = data.begin(), limit = data.end();
  Coord upperLimit, lowerLimit, last = result.back();

  for (++iter; iter != limit; ++iter) {
    // Add the new bounds
    Coord u = *iter, l = *iter, b = result.back();
    u.second += epsilon;
    l.second -= epsilon;

    // Check if we cut the error corridor
    if ((last != b) && ((cmpDevs(b, upperLimit, *iter) < 0) || (cmpDevs(b, lowerLimit, *iter) > 0))) {
      result.push_back(last);
      b = last;
    }

    // Update the error margins
    if ((last == b) || (cmpDevs(b, upperLimit, u) > 0))
      upperLimit = u;
    if ((last == b) || (cmpDevs(b, lowerLimit, l) < 0))
      lowerLimit = l;

    // And remember the current point
    last = *iter;
  }

  // Add the last point
  result.push_back(*(--data.end()));

  return result;
}

static std::vector<Coord> compressFunc(const std::vector<Coord> &func, unsigned desiredSize)
// Compress to the desired size
{
  // Relax a bit to speed up compression
  unsigned maxSize = desiredSize + (desiredSize / 100);
  unsigned minSize = desiredSize - (desiredSize / 100);

  // Fits?
  if (func.size() <= maxSize)
    return func;

  // No, binary search for a good enough epsilon (== middle) in (0, func.size())
  unsigned left = 0, right = func.size();
  while (left < right) {
    unsigned middle = (left + right) / 2;
    std::vector<Coord> candidate = tautString(func, func.back().second, middle);
    if (candidate.size() < minSize) {
      right = middle;
    } else if (candidate.size() > maxSize) {
      left = middle + 1;
    } else {
      return candidate;
    }
  }

  // Final call, this is the best we could get
  return tautString(func, func.back().second, left);
}

static std::vector<double> computeSlopes(const std::vector<Coord> &spline) __attribute__((unused));
static std::vector<double> computeSlopes(const std::vector<Coord> &spline)
// spare the divisions in the lookups by precomputing the slopes between each node of the spline
{
  uint32_t spline_size = spline.size();
  std::vector<double> slopes(spline_size - 1);
  for (uint32_t index = 0; index < spline_size - 1; ++index) {
    double dx = spline[index].first - spline[index + 1].first;
    double dy = spline[index].second - spline[index + 1].second;
    slopes[index] = dy / dx;
  }
  return slopes;
}

static std::vector<Coord> buildCdf(const std::vector<Geography::Point> &points) {
  std::vector<Coord> cdf;
  {
    unsigned pos = 0;
    double last = points.front().y;
    for (auto d : points) {
      if (d.y != last) {
        cdf.push_back({last, pos - 1});
        last = d.y;
      }
      pos++;
    }
    cdf.push_back({last, pos - 1});
  }
  return cdf;
}

struct CdfOnTheFlyInterfaceY {
  CdfOnTheFlyInterfaceY(const std::vector<Geography::Point> &points)
      : points(points) {
    assert(!points.empty());

    double curr = points[idx].y;
    while (idx + 1 < points.size() && curr == points[idx + 1].y) {
      idx++;
    }
    current = Coord{points[idx].y, idx};
  }

  const Coord &GetCurrent() {
    return current;
  }

  void Next() {
    assert(HasCurrent());

    idx++;
    if (idx == points.size()) {
      return;
    }

    double curr = points[idx].y;
    while (idx + 1 < points.size() && curr == points[idx + 1].y) {
      idx++;
    }
    current = Coord{points[idx].y, idx};
  }

  bool HasCurrent() {
    return idx < points.size();
  }

  uint32_t DataSize() {
    return points.size();
  }

 private:
  uint32_t idx = 0;
  Coord current;
  const std::vector<Geography::Point> &points;
};

// lowerLimit and upperLimit span a corridor, if the iter point is not in it we add it
// Otherwise, it is kind of on the same line and we dont need to add it
// On each iteration the corridor is narrowed
// Epsilon specifies the error bound: larger epsilon -> larger error corridor
static std::vector<Coord> tautString(CdfOnTheFlyInterfaceY &interface, double epsilon) {
  // Add the first point
  std::vector<Coord> result;
  result.reserve(interface.DataSize() * 2 / epsilon);
  if (!interface.HasCurrent())
    return result;
  result.push_back(interface.GetCurrent());
  interface.Next();
  if (!interface.HasCurrent())
    return result;

  // Taut the string
  Coord upperLimit, lowerLimit, last = result.back();

  while (interface.HasCurrent()) {
    // Add the new bounds
    Coord u = interface.GetCurrent(), l = interface.GetCurrent(), b = result.back();
    u.second += epsilon;
    l.second -= epsilon;

    // Check if we cut the error corridor
    if ((last != b) && ((cmpDevs(b, upperLimit, interface.GetCurrent()) < 0)
        || (cmpDevs(b, lowerLimit, interface.GetCurrent()) > 0))) {
      result.push_back(last);
      b = last;
    }

    // Update the error margins
    if ((last == b) || (cmpDevs(b, upperLimit, u) > 0))
      upperLimit = u;
    if ((last == b) || (cmpDevs(b, lowerLimit, l) < 0))
      lowerLimit = l;

    // And remember the current point
    last = interface.GetCurrent();
    interface.Next();
  }

  // Add the last point
  result.push_back(interface.GetCurrent());

  return result;
}
// ---------------------------------------------------------------------------------------------------------------------
}
