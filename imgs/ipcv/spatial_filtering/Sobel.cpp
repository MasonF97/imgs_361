/** Implementation file for image filtering
 *
 *  \file ipcv/spatial_filtering/Filter2D.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 20 Sep 2018
 */

#include "Sobel.h"

#include <iostream>
#include "imgs/ipcv/spatial_filtering/Filter2D.h"

using namespace std;

namespace ipcv {

/** Does a sobel gradient
 *
 *  \param[in] src          source cv::Mat of CV_8UC3
 *  \param[out] dst         destination cv::Mat of ddepth type
 *  \param[in] delta        optional value added to the filtered pixels
 *                          before storing them in dst
 *  \param[in] border_mode  pixel extrapolation method
 *  \param[in] border_value value to use for constant border mode
 */
bool Sobel(const cv::Mat& src, cv::Mat& dst,
                const int delta,
                const BorderMode border_mode,
                uint8_t border_value) {

  // Create the kernel to get the x gradient
  cv::Mat sobel_x_kernel = (cv::Mat_<float>(3, 3) <<
    -1, 0, 1,
    -2, 0, 2,
    -1, 0, 1);

  // Create the kernel to get the y gradient
  cv::Mat sobel_y_kernel = (cv::Mat_<float>(3, 3) <<
    -1, -2, -1,
    0,  0,  0,
    1,  2,  1);

  // Set the anchor to be in the center
  cv::Point anchor;
  anchor.x = -1;
  anchor.y = -1;
  cv::Mat grad_x, grad_y;
  // Use the Filter2D function to get the gradients
  ipcv::Filter2D(src, grad_x, CV_8UC3, sobel_x_kernel, anchor, delta, border_mode, border_value);
  ipcv::Filter2D(src, grad_y, CV_8UC3, sobel_y_kernel, anchor, delta, border_mode, border_value);

  // Output should only be one channel
  dst = cv::Mat::zeros(src.size(), CV_8UC1);
  for (int y = 0; y < src.rows; y++) {
    for (int x = 0; x < src.cols; x++) {
        // Use one of the channels
        uchar gx = grad_x.at<cv::Vec3b>(y, x)[0];
        uchar gy = grad_y.at<cv::Vec3b>(y, x)[0];

        // Get the absolute value
        int abs_val = std::abs((int)gx) + std::abs((int)gy);

        dst.at<uchar>(y, x) = cv::saturate_cast<uchar>(abs_val);
    }
  }
  return true;
}
}
