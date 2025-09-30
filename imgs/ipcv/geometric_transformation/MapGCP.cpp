/** Implementation file for finding source image coordinates for a source-to-map
 *  remapping using ground control points
 *
 *  \file ipcv/geometric_transformation/MapGCP.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 15 Sep 2018
 */

#include "MapGCP.h"

#include <iostream>
#include <vector>
#include <utility>

#include <Eigen/Dense>
#include <opencv2/core.hpp>

using namespace std;

namespace ipcv {

vector<std::pair<int, int>> get_polynomial_pairs(int order){
  vector<std::pair<int, int>> polynomial_pairs;
  for (int y = 0; y<= order; y++){
    for (int x = 0; x<= order; x++){
      if (x+y <= order){
        polynomial_pairs.push_back({x,y});
      }
    }
  }
  return polynomial_pairs;
}

/** Find the source coordinates (map1, map2) for a ground control point
 *  derived mapping polynomial transformation
 *
 *  \param[in] src   source cv::Mat of CV_8UC3
 *  \param[in] map   map (target) cv::Mat of CV_8UC3
 *  \param[in] src_points
 *                   vector of cv::Points representing the ground control
 *                   points from the source image
 *  \param[in] map_points
 *                   vector of cv::Points representing the ground control
 *                   points from the map image
 *  \param[in] order  mapping polynomial order
 *                      EXAMPLES:
 *                        order = 1
 *                          a0*x^0*y^0 + a1*x^1*y^0 +
 *                          a2*x^0*y^1
 *                        order = 2
 *                          a0*x^0*y^0 + a1*x^1*y^0 + a2*x^2*y^0 +
 *                          a3*x^0*y^1 + a4*x^1*y^1 +
 *                          a5*x^0*y^2
 *                        order = 3
 *                          a0*x^0*y^0 + a1*x^1*y^0 + a2*x^2*y^0 + a3*x^3*y^0 +
 *                          a4*x^0*y^1 + a5*x^1*y^1 + a6*x^2*y^1 +
 *                          a7*x^0*y^2 + a8*x^1*y^2 +
 *                          a9*x^0*y^3
 *  \param[out] map1  cv::Mat of CV_32FC1 (size of the destination map)
 *                    containing the horizontal (x) coordinates at which to
 *                    resample the source data
 *  \param[out] map2  cv::Mat of CV_32FC1 (size of the destination map)
 *                    containing the vertical (y) coordinates at which to
 *                    resample the source data
 */
bool MapGCP(const cv::Mat src, const cv::Mat map,
            const vector<cv::Point> src_points,
            const vector<cv::Point> map_points, const int order,
            cv::Mat& map1, cv::Mat& map2) {
  // Use a psuedo inverse for x
  // a = (Xt X)^-1 * Xt Y
  // X transpose * X gives a square matrix
  // Then get the inverse
  // Then get multiply it by the transpose again
  // Then multiply by the Y
  /*
  The sums of the exponents are not greater than 2, so it's a second order polynomial
  x' = a0 x0 y0 + a1 x1 y0 + a2 x0 y1 + a3 x1 y1 + a4 x2 y0 + a5 x0 y2
  Can make this using 2 loops
  x0 y0
  x0 y1 x1 y0
  x0 y2 x1  y1 x2 y0
  (Reverse the order)

  Model matrix = 6 terms in the equation, so 6 terms in the model equation
  X = 1 x1 y1 x1y1 x2 y2
  This is multiple linear least squares regression
  Observation matrix Y': Two matrices x' and y'

  Solve for A's and B's
  a = (Xt X)^-2 * Xt Yx'
  b = (Xt X)^-2 * Xt Yy'
  */

  // Construct model matrix

  // Get polynomials
  vector<std::pair<int, int>> polynomial_pairs = get_polynomial_pairs(order);

  for (std::pair<int, int> pair : polynomial_pairs){
    cout << pair.first << " " << pair.second << endl;
  }

  // Construct observation matrices



  return true;
}
}
