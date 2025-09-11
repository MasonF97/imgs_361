/** Implementation file for image enhancement using linear histogram statistics
 *
 *  \file ipcv/histogram_enhancement/LinearLut.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 3 Sep 2018
 */

#include "LinearLut.h"

#include <iostream>

#include "imgs/ipcv/utils/Utils.h"

using namespace std;

namespace ipcv {

// Clips a certain percent of the top and bottom of a cdf
void _clip_values_from_cdf(cv::Mat& cdf_row, int percentage) {
  // Calculate the upper and lower percentage thresholds
  float percentage_float = percentage / 100.0f;
  float lower_thresh = percentage_float / 2;
  float upper_thresh = 1 - (percentage_float / 2);

  // Go through each value
  for (int i = 0; i < 256; i++) {
    double& val = cdf_row.at<double>(i);
    // If the value is less than the lower threshold, set it to 0
    if (val <= lower_thresh) {
      val = 0.0;
    // If the value is greater than the upper threshold, set it to 1
    } else if (val >= upper_thresh) {
      val = 1.0;
    }
  }
}

// Maps lut values
void _create_lut_from_cdf(const cv::Mat& cdf_row, std::vector<uchar>& lut_values) {
  // Resize the lut_values so there isn't a segmentation fault
  lut_values.resize(256);
  // For each value, multiply it by 255 to get the value that should be in the LUT
  for (int i = 0; i < 256; ++i) {
    double val = cdf_row.at<double>(i);
    // use round and clamp to make sure the value is valid
    lut_values[i] = static_cast<uchar>(std::round(std::clamp(val * 255.0, 0.0, 255.0)));
  }
}

/** Create a 3-channel (color) LUT using linear histogram enhancement
 *
 *  \param[in] src          source cv::Mat of CV_8UC3
 *  \param[in] percentage   the total percentage to remove from the tails
 *                          of the histogram to find the extremes of the
 *                          linear enhancemnt function
 *  \param[out] lut         3-channel look up table in cv::Mat(3, 256)
 */
bool LinearLut(const cv::Mat& src, const int percentage, cv::Mat& lut) {

  // Create the histogram for the image
  cv::Mat histogram;
  ipcv::Histogram(src, histogram);
  // Create the cdf from the histogram
  cv::Mat cdf;
  ipcv::HistogramToCdf(histogram, cdf);

  // Create the empty output LUT
  lut = cv::Mat(3, 256, CV_8U);
  // Go through each row
  for (int c = 0; c < cdf.rows; ++c) {
    cv::Mat cdf_row = cdf.row(c);

    // If percentage is not 0, clip the values 
    if (percentage != 0){
      _clip_values_from_cdf(cdf_row, percentage);
    }
    
    // Get the values that will go into the lut
    std::vector<uchar> lut_values;
    _create_lut_from_cdf(cdf_row, lut_values);

    // Copy the lut_values into the lut
    std::memcpy(lut.ptr<uchar>(c), lut_values.data(), 256 * sizeof(uchar));
  }
  return true;
}
}
