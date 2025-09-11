/** Implementation file for computing a CDF from a histogram
 *
 *  \file ipcv/utils/HistogramToCdf.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 17 Mar 2018
 */

#include "HistogramToCdf.h"

#include "imgs/ipcv/utils/HistogramToPdf.h"

namespace ipcv {

void HistogramToCdf(const cv::Mat& h, cv::Mat& cdf) {
  cv::Mat pdf;
  HistogramToPdf(h, pdf);

  cdf.create(h.size(), CV_64F);

  for (int channel = 0; channel < pdf.rows; channel++) {
    double cumulative = 0.0;
    for (int digital_count_percent = 0; digital_count_percent < pdf.cols; digital_count_percent++) {
      // Update each value so that it represents the percentage up to that point
      cumulative += pdf.at<double>(channel, digital_count_percent);
      cdf.at<double>(channel, digital_count_percent) = cumulative;
    }
  }
}
}
