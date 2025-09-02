/** Implementation file for image quantization
 *
 *  \file ipcv/quantize/quantize.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 17 Mar 2018
 */

#include "Quantize.h"

#include <iostream>
#include <cmath>

using namespace std;


/** Determines whether a number is a power of 2
 *
 *  \param[in] number the number of being evaluated
 *
 *  \return a boolean indicating whether the number is a power of 1
 */
static int _is_power_of_two(const int number) {
  float log_2_double = log2(number);
  int log_2_int = static_cast<int>(log_2_double);
  return log_2_double == log_2_int;
}

/** Perform uniform grey-level quantization on a color image
 *
 *  \param[in] src                 source cv::Mat of CV_8UC3
 *  \param[in] quantization_levels the number of levels to which to quantize
 *                                 the image
 *  \param[out] dst                destination cv:Mat of CV_8UC3
 */
void Uniform(const cv::Mat& src, const int quantization_levels, cv::Mat& dst) {  

  // Might be more efficent to just always use division

  // Check whether the quantization level number is a power of two
  float quantization_denominator;
  if (_is_power_of_two(quantization_levels)){
    // Use bit shifting
    uint8_t quantization_levels_log_2_int = static_cast<int>(log2(quantization_levels));
    quantization_denominator = (256 >> quantization_levels_log_2_int);
  }
  else{
    // Use division
    // 256 must be a float so that a float is returned
    quantization_denominator = (256.0 / quantization_levels);
  }

  // Go th
  for (int r = 0; r < src.rows; r++) {
    for (int c = 0; c < src.cols; c++) {
      uint8_t src_value = src.at<cv::Vec3b>(r, c)[0];
      uint8_t quantized_value = src_value / quantization_denominator;
      dst.at<cv::Vec3b>(r, c) = cv::Vec3b(quantized_value, quantized_value, quantized_value);
    }
  }
}

/** Perform improved grey scale quantization on a color image
 *
 *  \param[in] src                 source cv::Mat of CV_8UC3
 *  \param[in] quantization_levels the number of levels to which to quantize
 *                                 the image
 *  \param[out] dst                destination cv:Mat of CV_8UC3
 */
void Igs(const cv::Mat& src, const int quantization_levels, cv::Mat& dst) {

  // Insert your code here

}

// create a copy of this repo to use
// fist get the env figured out
// run it as is, to make sure everyhtong works


// many functions
// use a base 2 logarithm to check if it's a power of two
// if it is a base 2 logarithmn use bit shifitng, otherwise use the method discussed in class
// a Mat of CV_8UC3 is an array of digital counts
// Only need to submit this Quantize.cpp file

// Look in ipcv/utils/ApplyLut.cpp


namespace ipcv {

bool Quantize(const cv::Mat& src, const int quantization_levels,
              const QuantizationType quantization_type, cv::Mat& dst) {
  dst.create(src.size(), src.type());

  switch (quantization_type) {
    case QuantizationType::uniform:
      Uniform(src, quantization_levels, dst);
      break;
    case QuantizationType::igs:
      Igs(src, quantization_levels, dst);
      break;
    default:
      cerr << "Specified quantization type is unsupported" << endl;
      return false;
  }

  return true;
}
}
