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

  // Insert your code here
  int dest_height = src.rows;
  int dest_width = src.cols;

  int kernel_anchor_row;
  int kernel_anchor_col;
  if (anchor.x == -1 && anchor.y == -1){
    kernel_anchor_row = (kernel.rows -1) / 2;
    kernel_anchor_col = (kernel.cols -1) / 2;
  }else{
    kernel_anchor_row = anchor.y;
    kernel_anchor_col = anchor.x;
  }
  // 256 -1  - (3 -1) - 2

  int max_x = (dest_width - 1) - ((kernel.cols -1) - kernel_anchor_col);
  int min_x = kernel_anchor_col;
  int max_y = (dest_height - 1) - ((kernel.rows -1) - kernel_anchor_row);
  int min_y = kernel_anchor_row;

  auto in_bounds = [max_x, min_x, max_y, min_y](int x, int y) {
    if (x < min_x || x > max_x){
      // cout << "y out of bounds" << endl;
      return false;
    }else if(y < min_y || y > max_y){
      // cout << "y out of bounds" << endl;
      return false;
    }else{
      return true;
    }
  };

  auto apply_kernel = [src, kernel, kernel_anchor_col, kernel_anchor_row](int x, int y){
    cv::Vec3f result(0,0,0);
    for (int s = (-1 * kernel_anchor_col); s < kernel_anchor_col + 1; s++){
      for (int t = (-1 * kernel_anchor_row); t < kernel_anchor_row +1; t++){
        // result += (kernel.at<float>(t + kernel_anchor_row, s + kernel_anchor_col) * src.at<cv::Vec3b>(y+t, x+s));
        cv::Vec3b pixel = src.at<cv::Vec3b>(y + t, x + s);
        float k = kernel.at<float>(t + kernel_anchor_row, s + kernel_anchor_col);
        result[0] += k * pixel[0];
        result[1] += k * pixel[1];
        result[2] += k * pixel[2];
      }
    }
    return cv::Vec3b(
      cv::saturate_cast<uchar>(result[0]),
      cv::saturate_cast<uchar>(result[1]),
      cv::saturate_cast<uchar>(result[2])
    );
  };

  cv::Vec3b delta_vector(delta, delta, delta);
  
  dst = cv::Mat::zeros(dest_height, dest_width, ddepth);

  for (int y = 0; y< dest_height; y++){
    for (int x = 0; x< dest_width; x++){
      if (in_bounds(x, y)){
        dst.at<cv::Vec3b>(y, x) = apply_kernel(x, y) + delta_vector;
      }else{
        // update for multiple depths
        // if(border_mode == BorderMode::CONSTANT)
        dst.at<cv::Vec3b>(y, x) = cv::Vec3b (border_value, border_value, border_value) + delta_vector;
      }

    }
  }
  return true;
}
}
