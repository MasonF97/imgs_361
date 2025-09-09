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

std::vector<uchar> _create_lut_for_color(const cv::Mat& channel, int percentage){
  cv::Mat flattened_channel = channel.reshape(1,1);
  std::vector<int> histogram(256, 0); 
  float total_pixels = static_cast<float>(flattened_channel.total()); // total number of pixels
  const uint8_t* data_ptr = flattened_channel.ptr<uint8_t>(0);
  for (size_t i = 0; i < total_pixels; ++i) {
    histogram[data_ptr[i]]++;
  }

  std::vector<int> cdf(256, 0);
  cdf[0] = histogram[0];
  for (int i = 1; i < 256; ++i) {
      cdf[i] = cdf[i - 1] + histogram[i];
  }

  // for (const auto& num : cdf) {
  //   std::cout << num << " ";
  // }

  // Step 2: Normalize the CDF to [0, 1] (or [0, 255] if needed)
  
  std::cout << total_pixels << std::endl;
  std::vector<float> cdf_normalized(256, 0.0f);

  for (int i = 0; i < 256; ++i) {
      cdf_normalized[i] = static_cast<float>(cdf[i]) / total_pixels;
  }

  std::cout << "Last CDF value: " << cdf[255] << std::endl;
  std::cout << "Total pixels  : " << total_pixels << std::endl;

  // Find the clip values
  float percentage_float = percentage / 100.0f;
  float lower_thresh = percentage_float / 2;
  float upper_thresh = 1 - (percentage_float / 2);

  // Set high values to 255, set low values to 0

  for (int i = 0; i< 256;i++){
    if (cdf_normalized[i] <= lower_thresh){
      // std::cout << "Found lower_clip_index at: " << i << std::endl;
      cdf_normalized[i] = 0;
    }
    else if (cdf_normalized[i] >= upper_thresh){
      cdf_normalized[i] = 1;
    }
  }

  // std::cout << "123" << std::endl;
  
  std::vector<uchar> lut_values(256, 0);
  // Multiply by 255 to get LUT value
  for (int i = 0; i< 256;i++){
    lut_values[i] = static_cast<uchar>(cdf_normalized[i] * 255);
  }
  // std::cout << "1234" << std::endl;

  // for (const auto& num : lut_values) {
  //   std::cout << num << " ";
  // }
  // std::cout << "12345" << std::endl;

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
  // To run: bin/histogram_enhance -i input_image -o destination_filename



  std::vector<cv::Mat> channels;
  cv::split(src, channels);
  // std::cout << channels.size() << std::endl;
  // std::cout << channels[0].size() << std::endl;

  // Create the empty output LUT
  lut = cv::Mat(3, 256, CV_8U);

  // figure out how to set the lut to the returned vectors
  for (uint8_t channel_index = 0;channel_index < channels.size(); channel_index++){
    std::vector<uchar> channel_lut_vector = _create_lut_for_color(channels[channel_index], percentage);
    // cv::Mat(channel_lut_vector, true).copyTo(lut.row(channel_index));
    // lut.row(channel_index) = cv::Mat(channel_lut_vector, true).clone(); 
    std::memcpy(lut.ptr(channel_index), channel_lut_vector.data(), 256);
  }
  // Print LUT
  // for (int row = 0; row < lut.rows; ++row) {
  //     std::cout << "Row " << row << ": [";
  //     for (int col = 0; col < lut.cols; ++col) {
  //         std::cout << static_cast<int>(lut.at<uchar>(row, col));
  //         if (col != lut.cols - 1) std::cout << ", ";
  //     }
  //     std::cout << "]" << std::endl;
  // }
  // Watch the lectures



  return true;
}
}
