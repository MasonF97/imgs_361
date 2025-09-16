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

  /*
  MATH
  We only need to compute thresholds where 0 < omega(k) < 1

  k is a index, 0-255

  omega(k) = sum from 0 to k {p(i)}

  mu(k) = sum from 0 to k {ip(i)}

  omega(c0) = omega(k)
  omega(c1) = 1 - omega(k)

  muT = sum from 0 to 255 {ip(i)}

  mu(c0) = mu(k)/omega(k)
  mu(c1) = (mu(T) - mu(k)) / (1 - omega(k))

  variance(c0) = sum from 0 to k {((i - mu(c0))^2)p(i)/omega(c0)}
  variance(c1) = sum from k+1 to 255 {((i - mu(c1))^2)p(i)/omega(c1)}

  variance(W) = (omega(c0)*variance(c0)) + (omega(c1)*variance(c1))
  variance(B) = omega(c0)*omega(c1)*((mu(c1) - mu(c0))^2)
  variance(T) = sum from 0 to 255 {((i - mu(T))^2)*p(i)}

  lambda = variance(B) / variance(W)
  kappa = variance(T) / variance(W)
  nu = variance(B) / variance(T)

  kappa = lamda + 1
  nu = lambda / (lambda + 1)

  nu(k) = variance(B)(k) / variance(T)
  variance(B)(k) = (((mu(T) * omega(k)) - mu(k))^2) / (omega(k)*(1 - omega(k)))

  The optimal threshold k* is the maximum value of variance(B)(k)
  */

int _get_threshold_for_color(cv::Mat& pdf){
  int length_of_pdf = 256;
  double current_max_inter_class_variance = 0.0;
  int current_max_k = 0;

  double omega_k = 0.0;
  double mu_k = 0.0;
  // total mean
  double mu_T = 0.0;
  for (int i = 0; i<length_of_pdf;i++){
    mu_T += i*pdf.at<double>(0,i);
  }

  for (int k = 0; k < length_of_pdf; k++) {
    // Add pdf(i) to omega_k during each iteration
    omega_k += pdf.at<double>(0, k);
    
    // if omega is not greater than 0 or is not smaller than 1, we can ignore it
    if (!(omega_k > 0 && omega_k < 1)){
      continue;
    }

    // mean of k
    mu_k += k*pdf.at<double>(0,k);

    double inter_class_variance = (pow(((mu_T * omega_k) - mu_k), 2)) / (omega_k * (1 - omega_k));

    if (inter_class_variance > current_max_inter_class_variance){
      current_max_inter_class_variance = inter_class_variance;
      current_max_k = k;
    }
  }
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

  // convert to histogram
  cv::Mat histogram;
  Histogram(src, histogram);

  // convert to PDF
  cv::Mat pdf;
  HistogramToPdf(histogram, pdf);

  for(int i = 0; i <3; i++){
    cv::Mat single_color_pdf = pdf.row(i).clone();
    threshold[i] = static_cast<uchar>(_get_threshold_for_color(single_color_pdf));
  }
  return true;
}
}
