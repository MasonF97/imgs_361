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
  if (radius < 0) {
    actual_radius = static_cast<int>(2 * sigma_distance);
  } else {
    actual_radius = radius;
  }

  // Precompute the spatial weights because they'll always be the same
  cv::Mat spatialKernel(2 * actual_radius + 1, 2 * actual_radius + 1, CV_32F);
  for (int i = -actual_radius; i <= actual_radius; i++) {
    for (int j = -actual_radius; j <= actual_radius; j++) {
      float dist2 = static_cast<float>(i * i + j * j);
      spatialKernel.at<float>(i + actual_radius, j + actual_radius) = 
          std::exp(-dist2 / (2 * sigma_distance * sigma_distance));
    }
  }

  // create a bool to indicate whether or not it's a color image
  // Check if all channels are identical to determine whether or not it's greyscale
  std::vector<cv::Mat> chans;
  cv::split(src, chans);
  bool is_greyscale = (cv::countNonZero(chans[0] != chans[1]) == 0 && cv::countNonZero(chans[0] != chans[2]) == 0);
  
  cv::Mat lab_or_gray;
  if (is_greyscale) {
    // if it's greyscale, just use src
    lab_or_gray = src.clone();  
  } else {
    // If it's a color image, convert to LAB
    cv::cvtColor(src, lab_or_gray, cv::COLOR_BGR2Lab);
  }
  
  // Create an array of channels that are floats
  std::vector<cv::Mat> float_channels;
  // Split up the channels and convert to float
  cv::split(lab_or_gray, float_channels);
  for (auto& ch : float_channels)
    ch.convertTo(ch, CV_32F);

  // Either use the L channel (or single grayscale channel) to compute the weights
  cv::Mat filter_channel = float_channels[0];
  
  // Create channels to store the float outputs
  // These will eventually be merged to create dst
  std::vector<cv::Mat> dst_channels(float_channels.size());
  for (size_t c = 0; c < float_channels.size(); c++) {
    dst_channels[c] = cv::Mat::zeros(src.size(), CV_32F);
  }

  // Go through each pixel
  for (int y = 0; y < rows; y++) {
    for (int x = 0; x < cols; x++) {
      float center = filter_channel.at<float>(y, x);
      float weight_sum = 0.0f;
      std::vector<float> pixel_sums(float_channels.size(), 0.0f);

      // Iterate through the filter
      for (int i = -actual_radius; i <= actual_radius; i++) {
        for (int j = -actual_radius; j <= actual_radius; j++) {
          int ny = y + i;
          int nx = x + j;

          float neighbor;
          // Check if the pixel is out of bounds
          if (ny < 0 || ny >= rows || nx < 0 || nx >= cols) {
            if (border_mode == BorderMode::CONSTANT) {
              neighbor = border_value;
            } else {
              // replicate mode
              int clamped_y = std::clamp(ny, 0, rows - 1);
              int clamped_x = std::clamp(nx, 0, cols - 1);
              neighbor = filter_channel.at<float>(clamped_y, clamped_x);
            }
          } else {
            neighbor = filter_channel.at<float>(ny, nx);
          }

          // Compute the range weight
          float range_diff = neighbor - center;
          float range_weight = std::exp(-(range_diff * range_diff) / (2 * sigma_range * sigma_range));

          // Compute the spatial weight and combine it with the range weight
          float spatial_weight = spatialKernel.at<float>(i + actual_radius, j + actual_radius);
          float weight = spatial_weight * range_weight;

          // Apply to each channel
          for (size_t c = 0; c < float_channels.size(); c++) {
            float channel_neighbor;
            // check if the neighbor is out of bounds
            if (ny < 0 || ny >= rows || nx < 0 || nx >= cols) {
              if (border_mode == BorderMode::CONSTANT) {
                channel_neighbor = border_value;
              } else {
                int clamped_y = std::clamp(ny, 0, rows - 1);
                int clamped_x = std::clamp(nx, 0, cols - 1);
                channel_neighbor = float_channels[c].at<float>(clamped_y, clamped_x);
              }
            } else {
              channel_neighbor = float_channels[c].at<float>(ny, nx);
            }
            // multiply the neighbor by the weight and then add it to the pixel sum
            pixel_sums[c] += channel_neighbor * weight;
          }
          weight_sum += weight;
        }
      }
      
      // Normalize and store the results for all channels
      for (size_t c = 0; c < float_channels.size(); c++) {
        dst_channels[c].at<float>(y, x) = pixel_sums[c] / weight_sum;
      }
    }
  }

  // Convert back to int
  for (auto& ch : dst_channels){
    ch.convertTo(ch, CV_8U);
  }
  // merge the output channels into dst
  cv::merge(dst_channels, dst);

  // Convert back to BGR if it's a color image
  if (!is_greyscale) {
    cv::cvtColor(dst, dst, cv::COLOR_Lab2BGR);
  }

  return true;
}
}