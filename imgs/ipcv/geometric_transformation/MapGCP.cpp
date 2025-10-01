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
  // Generate all x, y exponent pairs for the polynomial order
 // Each pair corresponds to a term x^i * y^j
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

double evaluate_polynomial(vector<pair<int, int>> polynomial_pairs, Eigen::MatrixXd coeff_matrix, double x, double y){
  // Evaluate the polynomial at x, y using the coefficients in the coeff_matrix 
  double sum = 0;
  int size_of_polynomial_pairs = static_cast<int>(polynomial_pairs.size());
  for( int j = 0; j < size_of_polynomial_pairs; j++){
    sum += coeff_matrix(j,0)*(pow(x, polynomial_pairs[j].first) * pow(y, polynomial_pairs[j].second));
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
  // Get the polynomial term exponents
  vector<std::pair<int, int>> polynomial_pairs = get_polynomial_pairs(order);

  // Construct the model matrix
  int model_matrix_rows = static_cast<int>(src_points.size());
  int model_matrix_cols = static_cast<int>(polynomial_pairs.size());
  Eigen::MatrixXd model_matrix(model_matrix_rows, model_matrix_cols);
  for (int i = 0; i < model_matrix_rows; i++){
    for( int j = 0; j < model_matrix_cols; j++){
      // calculate the value using the exponents in polynomial_pairs
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

  // Solve for the a (x) and b (y) coefficient matrices using the least squares solution
  // Uses the psuedo inverse
  Eigen::MatrixXd a_coefficient_matrix = (model_matrix.transpose() * model_matrix).inverse() * model_matrix.transpose() * x_observation_matrix;
  Eigen::MatrixXd b_coefficient_matrix = (model_matrix.transpose() * model_matrix).inverse() * model_matrix.transpose() * y_observation_matrix;

  // Create the maps using the map image's size
  int dest_height = map.rows;
  int dest_width = map.cols;

  map1 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  map2 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);

  for(int y = 0; y < dest_height; y++){
    for (int x = 0; x < dest_width; x++){
      // Get the x source coordinate 
      map1.at<float>(y, x) = evaluate_polynomial(polynomial_pairs, a_coefficient_matrix, x, y);
      // Get the x source coordinate 
      map2.at<float>(y, x) = evaluate_polynomial(polynomial_pairs, b_coefficient_matrix, x, y);
    }
  }

  return true;
}
}
