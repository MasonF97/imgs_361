/** Implementation file for mapping a source quad on to a target quad
 *
 *  \file ipcv/geometric_transformation/MapPolar.cpp
 */

#include "MapPolar.h"

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <opencv2/core.hpp>

using namespace std;

namespace ipcv {

/** Find the source coordinates (map1, map2) for a polar mapping
 *
 *  \param[in] src       source cv::Mat of CV_8UC3
 *  \param[in] use_log   whether or not to do a log polar mapping
 *  \param[out] map1     cv::Mat of CV_32FC1 (size of the destination map)
 *                       containing the horizontal (x) coordinates at
 *                       which to resample the source data
 *  \param[out] map2     cv::Mat of CV_32FC1 (size of the destination map)
 *                       containing the vertical (y) coordinates at
 *                       which to resample the source data
 */
bool MapPolar(const cv::Mat src, const bool use_log, cv::Mat& map1,
            cv::Mat& map2) {

  int src_width = src.cols;
  int src_height = src.rows;

  // Find the center coords
  double cx = src_width / 2.0;
  double cy = src_height / 2.0;


  // get destination dimensions
  int dest_width = src.cols;
  int dest_height = src.rows;

  // Find the log base used to find r if doing a log polar mapping
  double log_base = 1.0;
  if (use_log) {
    log_base = std::pow(10, (std::log10(dest_width)) / static_cast<double>(dest_height));
  }

  map1 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  map2 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);

  for (int col_idx = 0; col_idx < dest_width; ++col_idx) {
    // Find the radius
    double r = use_log ? std::pow(log_base, col_idx) - 1 : static_cast<double>(col_idx);
    // double r = use_log ? std::pow(log_base, col_idx) - 1 : (static_cast<double>(col_idx) / (dest_width - 1)) * max_radius;
    
    for (int row_idx = 0; row_idx < dest_height; ++row_idx) {
      // theta, going CW, starting in the north, bottom to top
      double theta = 2.0 * M_PI * (static_cast<double>(row_idx) / (dest_height));
      // Find the source coord
      double src_x = cx + r * std::cos(theta);
      double src_y = cy + r * std::sin(theta);

      // Set the source coord
      map1.at<float>((dest_height - 1) - row_idx, col_idx) = static_cast<float>(src_x);
      map2.at<float>((dest_height - 1) - row_idx, col_idx) = static_cast<float>(src_y);
    }
  }

  return true;
}
}
