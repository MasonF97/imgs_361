/** Implementation file for finding Otsu's threshold
 *
 *  \file ipcv/otsus_threshold/OtsusThreshold.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 8 Sep 2018
 */

#include "OtsusThreshold.h"

#include <iostream>

#include <cmath>

#include "imgs/ipcv/utils/Utils.h"

using namespace std;

namespace ipcv {

// Get the optimal threshold value for a color using it's pdf and Otsu's method
int _get_threshold_for_color(cv::Mat& pdf){
  int length_of_pdf = 256;
  // Track the current k with the highest inter class variance
  double current_max_inter_class_variance = 0.0;
  int current_max_k = 0;

  double omega_k = 0.0;
  double mu_k = 0.0;
  // calculate the total mean
  double mu_T = 0.0;
  for (int i = 0; i<length_of_pdf;i++){
    mu_T += i*pdf.at<double>(0,i);
  }
  // calculate the inter class variance of each possible value of k
  for (int k = 0; k < length_of_pdf; k++) {
    // Add pdf(i) to omega_k during each iteration
    omega_k += pdf.at<double>(0, k);
    
    // if omega is not greater than 0 or is not smaller than 1, we can ignore it
    if (!(omega_k > 0 && omega_k < 1)){
      continue;
    }

    // mean of k
    mu_k += k*pdf.at<double>(0,k);

    // Calculate the inter class variance using the formula in Otsu's paper
    double inter_class_variance = (pow(((mu_T * omega_k) - mu_k), 2)) / (omega_k * (1 - omega_k));

    // Check if the variance is greater than the current highest variance.
    if (inter_class_variance > current_max_inter_class_variance){
      current_max_inter_class_variance = inter_class_variance;
      current_max_k = k;
    }
  }
  // The k-value with the highest inter class variance will be the optimal threshold according to Otsu's paper
  return current_max_k;
}

/** Find Otsu's threshold for each channel of a 3-channel (color) image
 *
 *  \param[in] src          source cv::Mat of CV_8UC3
 *  \param[out] threshold   threshold values for each channel of a 3-channel
 *                          color image in cv::Vec3b
 */
bool OtsusThreshold(const cv::Mat& src, cv::Vec3b& threshold) {
  threshold = cv::Vec3b();

  // convert src to histogram
  cv::Mat histogram;
  Histogram(src, histogram);

  // convert the histogram to a PDF
  cv::Mat pdf;
  HistogramToPdf(histogram, pdf);

  // For each color, calculate the threshold and add it to the thresholds vector
  for(int i = 0; i <3; i++){
    cv::Mat single_color_pdf = pdf.row(i).clone();
    threshold[i] = static_cast<uchar>(_get_threshold_for_color(single_color_pdf));
  }
  return true;
}
}
