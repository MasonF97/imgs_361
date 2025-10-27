/** Implementation file for bilateral filtering
 *
 *  \file ipcv/bilateral_filtering/BilateralFilter.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 29 Sep 2018
 */

#include "BilateralFilter.h"

#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;

namespace ipcv {

/** Bilateral filter an image
 *
 *  \param[in] src             source cv::Mat of CV_8UC3
 *  \param[out] dst            destination cv::Mat of ddepth type
 *  \param[in] sigma_distance  standard deviation of distance/closeness filter
 *  \param[in] sigma_range     standard deviation of range/similarity filter
 *  \param[in] radius          radius of the bilateral filter (if negative, use
 *                             twice the standard deviation of the distance/
 *                             closeness filter)
 *  \param[in] border_mode     pixel extrapolation method
 *  \param[in] border_value    value to use for constant border mode
 */
bool BilateralFilter(const cv::Mat& src, cv::Mat& dst,
                     const double sigma_distance, const double sigma_range,
                     const int radius, const BorderMode border_mode,
                     uint8_t border_value) {
  int rows = src.rows;
  int cols = src.cols;
  dst.create(src.size(), src.type());

  // Check if the radius is negative
  int actual_radius;
  if ( radius < 0){
    actual_radius = static_cast<int>(2 * sigma_distance);
  }else{
    actual_radius = radius;
  }

  // Precompute spatial Gaussian weights
  cv::Mat spatialKernel(2 * actual_radius + 1, 2 * actual_radius + 1, CV_32F);
  for (int i = -actual_radius; i <= actual_radius; i++) {
    for (int j = -actual_radius; j <= actual_radius; j++) {
      float dist2 = static_cast<float>(i * i + j * j);
      spatialKernel.at<float>(i + actual_radius, j + actual_radius) = std::exp(-dist2 / (2 * sigma_distance * sigma_distance));
    }
  }

  int num_of_channels = src.channels();

  cv::Mat src_LAB;
  std::vector<cv::Mat> channels;
  std::vector<cv::Mat> float_channels;
  bool color_img = false;
  if (num_of_channels == 3){
    // convert the src from BGR to LAB
    cv::cvtColor(src, src_LAB, cv::COLOR_BGR2Lab);
    cv::split(src_LAB, channels);
    color_img = true;
  }else{
    cv::split(src, channels);
  }
  cv::Mat filter_channel = channels[0];
  filter_channel.convertTo(filter_channel, CV_32F);
  cv::Mat dst_channel = cv::Mat::zeros(src.size(), CV_32F);

  // Go through each pixel
  for (int y = 0; y< rows; y++){
    for(int x = 0; x< cols; x++){
      float center = filter_channel.at<float>(y,x);
      float weight_sum = 0.0f;
      float pixel_sum = 0.0f;

      // Iterate through the filter
      for (int i = -actual_radius; i <= actual_radius; i++) {
        for (int j = -actual_radius; j <= actual_radius; j++) {
          int ny = y + i;
          int nx = x + j;

          float neighbor;
          // Check if the pixel is out of bounds
          if (ny < 0 || ny >= rows || nx < 0 || nx >= cols){
            if (border_mode == BorderMode::CONSTANT){
              neighbor = border_value;
            }else{
              // replicate mode
              int clamped_y = std::clamp(ny, 0, rows - 1);
              int clamped_x = std::clamp(nx, 0, cols - 1);
              neighbor = filter_channel.at<float>(clamped_y, clamped_x);
            }
          }else{
            neighbor = filter_channel.at<float>(ny, nx);
          }


          // Compute range weight
          float range_diff = (neighbor - center) / 100.0f;
          float range_weight = std::exp(-(range_diff * range_diff) / (2 * sigma_range * sigma_range));

          // Combine with spatial weight
          float spatial_weight = spatialKernel.at<float>(i + actual_radius, j + actual_radius);
          float weight = spatial_weight * range_weight;

          // Accumulate
          pixel_sum += neighbor * weight;
          weight_sum += weight;
        }
      }
      dst_channel.at<float>(y, x) = pixel_sum / weight_sum;

    }
  }
  dst_channel.convertTo(channels[0], CV_8U);
  cv::merge(channels, dst);

  if (color_img){
    cv::cvtColor(dst, dst, cv::COLOR_Lab2BGR);
  }

  return true;
}
}
