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

double get_prime_coord(vector<pair<int, int>> polynomial_pairs, Eigen::MatrixXd a_or_b, double x, double y){
  double sum = 0;
  int size_of_polynomial_pairs = static_cast<int>(polynomial_pairs.size());
  for( int j = 0; j < size_of_polynomial_pairs; j++){
    sum += a_or_b(j,0)*(pow(x, polynomial_pairs[j].first) * pow(y, polynomial_pairs[j].second));
  }
  return sum;
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

  // Get polynomials
  vector<std::pair<int, int>> polynomial_pairs = get_polynomial_pairs(order);

  // for (std::pair<int, int> pair : polynomial_pairs){
  //   cout << pair.first << " " << pair.second << endl;
  // }

  // Construct model matrix
  int model_matrix_rows = static_cast<int>(src_points.size());
  int model_matrix_cols = static_cast<int>(polynomial_pairs.size());
  // Does it need to be a double?
  // vector<vector<double>> model_matrix (model_matrix_rows, std::vector<double>(model_matrix_cols, 0));
  Eigen::MatrixXd model_matrix(model_matrix_rows, model_matrix_cols);
  for (int i = 0; i < model_matrix_rows; i++){
    for( int j = 0; j < model_matrix_cols; j++){
      model_matrix(i, j) = pow(map_points[i].x, polynomial_pairs[j].first) * pow(map_points[i].y, polynomial_pairs[j].second);
    }
  }

  // Construct observation matrices
  // Create x observation matrix
  Eigen::MatrixXd x_observation_matrix (model_matrix_rows, 1);
  // Create y observation matrix
  Eigen::MatrixXd y_observation_matrix (model_matrix_rows, 1);
  // Populate both matrices
  for (int i = 0; i< model_matrix_rows;i++){
    x_observation_matrix(i,0) = src_points[i].x;
    y_observation_matrix(i,0) = src_points[i].y;
  }

  Eigen::MatrixXd a_matrix = (model_matrix.transpose() * model_matrix).inverse() * model_matrix.transpose() * x_observation_matrix;
  Eigen::MatrixXd b_matrix = (model_matrix.transpose() * model_matrix).inverse() * model_matrix.transpose() * y_observation_matrix;


  int dest_height = map.rows;
  int dest_width = map.cols;

  map1 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  map2 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);

  for(int y = 0; y < dest_height; y++){
    for (int x = 0; x < dest_width; x++){
      // each polynomial + a or b
      map1.at<float>(y, x) = get_prime_coord(polynomial_pairs, a_matrix, x, y);
      map2.at<float>(y, x) = get_prime_coord(polynomial_pairs, b_matrix, x, y);
    }
  }

  return true;
}
}
