/** Implementation file for image filtering
 *
 *  \file ipcv/spatial_filtering/Filter2D.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 20 Sep 2018
 */

#include "Filter2D.h"

#include <iostream>

using namespace std;

namespace ipcv {

/** Correlates an image with the provided kernel
 *
 *  \param[in] src          source cv::Mat of CV_8UC3
 *  \param[out] dst         destination cv::Mat of ddepth type
 *  \param[in] ddepth       desired depth of the destination image
 *  \param[in] kernel       convolution kernel (or rather a correlation
 *                          kernel), a single-channel floating point matrix
 *  \param[in] anchor       anchor of the kernel that indicates the relative
 *                          position of a filtered point within the kernel;
 *                          the anchor should lie within the kernel; default
 *                          value (-1,-1) means that the anchor is at the
 *                          kernel center
 *  \param[in] delta        optional value added to the filtered pixels
 *                          before storing them in dst
 *  \param[in] border_mode  pixel extrapolation method
 *  \param[in] border_value value to use for constant border mode
 */
bool Filter2D(const cv::Mat& src, cv::Mat& dst, const int ddepth,
              const cv::Mat& kernel, const cv::Point anchor, const int delta,
              const BorderMode border_mode, const uint8_t border_value) {
  int dest_height = src.rows;
  int dest_width = src.cols;

  int kernel_anchor_row;
  int kernel_anchor_col;
  // Find the row and column of the anchor
  if (anchor.x == -1 && anchor.y == -1){
    kernel_anchor_row = (kernel.rows -1) / 2;
    kernel_anchor_col = (kernel.cols -1) / 2;
  }else{
    kernel_anchor_row = anchor.y;
    kernel_anchor_col = anchor.x;
  }
  // Find the min and max y, used for the isolated border mode
  int max_x = (dest_width - 1) - ((kernel.cols -1) - kernel_anchor_col);
  int min_x = kernel_anchor_col;
  int max_y = (dest_height - 1) - ((kernel.rows -1) - kernel_anchor_row);
  int min_y = kernel_anchor_row;
  // Create the destination Mat
  dst = cv::Mat::zeros(dest_height, dest_width, src.type());

  // Iterate through each destination pixel
  for (int y = 0; y< dest_height; y++){
    for (int x = 0; x< dest_width; x++){
      // Create a vector to store the results of adding up each square in the kernel
      cv::Vec3f result(0,0,0);
      // Check if the border mode is isolated, and if it is, check if the x,y are out of bounds
      if (border_mode == BorderMode::ISOLATED && (x < min_x || x > max_x || y < min_y || y > max_y)){
        dst.at<cv::Vec3b>(y, x) = cv::Vec3b(
        cv::saturate_cast<uchar>(border_value + delta),
        cv::saturate_cast<uchar>(border_value + delta),
        cv::saturate_cast<uchar>(border_value + delta));
        continue;  // Skip filtering
      }
      // Iterate through the kernel
      for (int s = 0; s < kernel.rows; s++){
        for (int t = 0; t < kernel.cols; t++){
          // Find the offset needed using the anchor position
          int y_offset = s - kernel_anchor_row;
          int x_offset = t - kernel_anchor_col;

          int src_y = y + y_offset;
          int src_x = x + x_offset;
          cv::Vec3b pixel;
          // Check whether the pixel is within the bounds of the image
          if (src_x >= 0 && src_x < dest_width && src_y >= 0 && src_y < dest_height){
            // If it is, get the pixel from the src image
            pixel = src.at<cv::Vec3b>(src_y, src_x);
          }else{
            // If it isn't check the border mdoe
            if(border_mode == BorderMode::CONSTANT){
              pixel = cv::Vec3b (border_value, border_value, border_value);
            }
            else if (border_mode == BorderMode::REPLICATE){
              // If border mode is replicate, clamp the coords to get the nearest pixel in the image
              int clamped_y = std::clamp(src_y, 0, dest_height - 1);
              int clamped_x = std::clamp(src_x, 0, dest_width - 1);
              pixel = src.at<cv::Vec3b>(clamped_y, clamped_x);
            }else{
              // Wrap mode
              // Add the dest_height and the modulo by the dest height, which will give us the correct wrap value              
              int wrapped_y = (src_y + dest_height) % dest_height;
              int wrapped_x = (src_x + dest_width) % dest_width;
      
              pixel = src.at<cv::Vec3b>(wrapped_y, wrapped_x);
            }
          }
          // Get the 'weight' from the kernel
          float k = kernel.at<float>(s, t);
          // Multiply the kernel value and the pixel value and add it to the result
          result[0] += k * pixel[0];
          result[1] += k * pixel[1];
          result[2] += k * pixel[2];
        }
      }
      // Set the destination pixel to the result pixel and add the delta value
      dst.at<cv::Vec3b>(y, x) = cv::Vec3b(
        cv::saturate_cast<uchar>(result[0] + delta),
        cv::saturate_cast<uchar>(result[1] + delta),
        cv::saturate_cast<uchar>(result[2] + delta)
      );

    }
  }
  return true;
}
}
