/** Implementation file for computing the PDF from a histogram
 *
 *  \file ipcv/utils/HistogramToPdf.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 17 Mar 2018
 */

#include "HistogramToPdf.h"

namespace ipcv {

void HistogramToPdf(const cv::Mat& h, cv::Mat& pdf) {
  pdf.create(h.size(), CV_64F);

  for (int row = 0; row < h.rows; row++) {
    // Get the total number of 'counts' in the row
    double row_sum = cv::sum(h.row(row))[0];
    for (int col = 0; col < h.cols; col++) {
      int count = h.at<int>(row, col);
      // For each value, divide it by the total number of values in the row, so that its between 0 and 1
      pdf.at<double>(row, col) = static_cast<double>(count) / row_sum;
    }
  }

}
}
