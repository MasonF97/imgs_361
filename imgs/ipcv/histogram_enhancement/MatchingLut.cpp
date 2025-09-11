/** Implementation file for image enhancement using histogram matching
 *
 *  \file ipcv/histogram_enhancement/MatchingLut.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 3 Sep 2018
 */

#include "MatchingLut.h"

#include <iostream>

#include "imgs/ipcv/utils/Utils.h"

using namespace std;

namespace ipcv {

// Return the index of the value closest to the target value in the row
int findClosestIndexMat(const cv::Mat& row, float target_value) {
  // Create a mat that has the difference between each value and the target value
  cv::Mat difference_mat;
  cv::absdiff(row, target_value, difference_mat);
  double min_val;
  cv::Point min_location;
  // Find the smallest value in the difference mat,
  // as the original value will the the 'closest' to the target value
  cv::minMaxLoc(difference_mat, &min_val, nullptr, &min_location, nullptr);
  // Return that index
  return min_location.x;
}

/** Create a 3-channel (color) LUT using histogram matching
 *
 *  \param[in] src   source cv::Mat of CV_8UC3
 *  \param[in] h     the histogram in cv:Mat(3, 256) that the
 *                   source is to be matched to
 *  \param[out] lut  3-channel look up table in cv::Mat(3, 256)
 */
bool MatchingLut(const cv::Mat& src, const cv::Mat& h, cv::Mat& lut) {
  // Calculate the histogram and cdf of the source image
  cv::Mat src_histogram;
  ipcv::Histogram(src, src_histogram);
  cv::Mat src_cdf;
  ipcv::HistogramToCdf(src_histogram, src_cdf);

  // Calculate the cdf of the target histogram
  cv::Mat target_cdf;
  ipcv::HistogramToCdf(h, target_cdf);

  // Create an empty LUT
  lut = cv::Mat(3, 256, CV_8U);

  for (int ch = 0; ch < 3; ++ch) {
    // get the target cdf row
    cv::Mat target_cdf_row = target_cdf.row(ch);
    for (int i = 0; i < 256; ++i) {
      // get the src cdf value at the index
      float src_val = static_cast<float>(src_cdf.at<double>(ch, i));
      // Find the index in the target_cdf row that's closest to the src_val
      int matchedIndex = findClosestIndexMat(target_cdf_row, src_val);
      // Set the value of the LUT at that index to the 'closest' index in the target cdf
      lut.at<uchar>(ch, i) = static_cast<uchar>(matchedIndex);
    }
  }
  return true;
}
}
