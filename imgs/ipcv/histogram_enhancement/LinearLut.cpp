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

void _create_histogram(const cv::Mat channel, int total_pixels, std::vector<int>& histogram){
  cv::Mat flattened_channel = channel.reshape(1,1);
  const uint8_t* data_ptr = flattened_channel.ptr<uint8_t>(0);
  for (int i = 0; i < total_pixels; i++) {
    histogram[data_ptr[i]]++;
  }
}


void _create_cdf(std::vector<int> histogram, std::vector<int>& cdf){
  cdf[0] = histogram[0];
  for (int i = 1; i < 256; i++) {
      cdf[i] = cdf[i - 1] + histogram[i];
  }
}

void _normalize_cdf(std::vector<int> cdf, int total_pixels, std::vector<float>& cdf_normalized){
  for (int i = 0; i < 256; i++) {
    cdf_normalized[i] = static_cast<float>(cdf[i]) / total_pixels;
  }
}

void _clip_values_from_cdf(std::vector<float>& cdf_normalized, int percentage){
  // Find the clip values
  float percentage_float = percentage / 100.0f;
  float lower_thresh = percentage_float / 2;
  float upper_thresh = 1 - (percentage_float / 2);

  // Set high values to 255, set low values to 0
  for (int i = 0; i< 256;i++){
    if (cdf_normalized[i] <= lower_thresh){
      cdf_normalized[i] = 0;
    }
    else if (cdf_normalized[i] >= upper_thresh){
      cdf_normalized[i] = 1;
    }
  }
}

void _create_lut_from_cdf(std::vector<float> cdf_normalized, std::vector<uchar>& lut_values){
  // Multiply by 255 to get LUT value
  for (int i = 0; i< 256;i++){
    lut_values[i] = static_cast<uchar>(cdf_normalized[i] * 255);
  }
}

std::vector<uchar> _create_lut_for_channel(const cv::Mat& channel, int percentage){
  // Find the total number of pixels
  int total_pixels = static_cast<int>(channel.total());

  // Create a histogram for the channel
  std::vector<int> histogram(256, 0);
  _create_histogram(channel, total_pixels, histogram); 

  // Create a cdf for the histogram
  std::vector<int> cdf(256, 0);
  _create_cdf(histogram, cdf);

  // Step 2: Normalize the CDF to [0, 1] (or [0, 255] if needed)
  // std::cout << total_pixels << std::endl;
  std::vector<float> cdf_normalized(256, 0.0f);
  _normalize_cdf(cdf, total_pixels, cdf_normalized);

  // Clip values
  if (percentage > 0) {
    _clip_values_from_cdf(cdf_normalized, percentage);
  }
  
  std::vector<uchar> lut_values(256, 0);
  _create_lut_from_cdf(cdf_normalized, lut_values);
  return lut_values;
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
  std::vector<cv::Mat> channels;
  cv::split(src, channels);
  // std::cout << channels.size() << std::endl;
  // std::cout << channels[0].size() << std::endl;

  // Create the empty output LUT
  lut = cv::Mat(3, 256, CV_8U);

  // figure out how to set the lut to the returned vectors
  for (uint8_t channel_index = 0;channel_index < channels.size(); channel_index++){
    std::vector<uchar> channel_lut_vector = _create_lut_for_channel(channels[channel_index], percentage);
    std::memcpy(lut.ptr(channel_index), channel_lut_vector.data(), 256);
  }
  // Print LUT
  for (int row = 0; row < lut.rows; ++row) {
      std::cout << "Row " << row << ": [";
      for (int col = 0; col < lut.cols; ++col) {
          std::cout << static_cast<int>(lut.at<uchar>(row, col));
          if (col != lut.cols - 1) std::cout << ", ";
      }
      std::cout << "]" << std::endl;
  }
  // Watch the lectures
  return true;
}
}
