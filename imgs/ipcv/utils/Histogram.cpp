/** Implementation file for computing an image histogram
 *
 *  \file ipcv/utils/Histogram.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 17 Mar 2018
 */

#include "Histogram.h"

namespace ipcv {

void Histogram(const cv::Mat& src, cv::Mat& h) {
  h = cv::Mat_<int>::zeros(3, 256);
  // Iterate through each pixel
  for (int y = 0; y < src.rows; y++) {
    for (int x = 0; x < src.cols; x++) {
      // Get the pixel with the digital count for each color
      cv::Vec3b digital_count_vector = src.at<cv::Vec3b>(y, x);
      // Add one to the histogram value/count at the corresponding index
      for (int c = 0; c < 3; c++) {
        int digital_count = digital_count_vector[c];
        h.at<int>(c, digital_count)++;
      }
    }
  }
}
}
