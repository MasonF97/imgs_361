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

/** Applies a transform lambda to each digital count in an image
 *
 *  \param[in] src         source cv::Mat of CV_8UC3
 *  \param[out] dst        destination cv:Mat of CV_8UC3
 *  \param[in] transform   a lambda that contains the calculations to 
 *                         perfrom uniform or igs quantization on a digital count
 */
static void _transform_all_digital_counts(const cv::Mat& src, cv::Mat& dst, std::function<uint8_t(uint8_t)> transform) {
  // Go through each digital count in the image
  for (int r = 0; r < src.rows; r++) {
    for (int c = 0; c < src.cols; c++) {
      uint8_t src_digital_count = src.at<cv::Vec3b>(r, c)[0];
      // Apply the transformation lambda to the digital count
      uint8_t quantized_digital_count = transform(src_digital_count);
      dst.at<cv::Vec3b>(r, c) = cv::Vec3b(quantized_digital_count, quantized_digital_count, quantized_digital_count);
    }
  }
}

/** Perform uniform grey-level quantization on a color image
 *
 *  \param[in] src                 source cv::Mat of CV_8UC3
 *  \param[in] quantization_levels the number of levels to which to quantize
 *                                 the image
 *  \param[out] dst                destination cv:Mat of CV_8UC3
 */
void Uniform(const cv::Mat& src, const int quantization_levels, cv::Mat& dst) {
  std::function<uint8_t(uint8_t)> uniform_transform;
  // Must be a float to avoid errors with odd quantization levels
  float quantization_denominator = (256.0 / quantization_levels);
  // Check whether the quantization level number is a power of two
  if (_is_power_of_two(quantization_levels)){
    // Determine the correct amount of bits to shift by
    uint8_t bit_shift = static_cast<uint8_t>(log2(quantization_denominator));
    // Create a lambda that uses bit shifting to perform uniform quantization
    uniform_transform = [bit_shift](uint8_t src_digital_count) {
      return src_digital_count >> bit_shift;
    };
  }
  else{
    // Create a lambda that uses division to perform uniform quantization
    uniform_transform = [quantization_denominator](uint8_t src_digital_count) {
      return static_cast<uint8_t>(src_digital_count / quantization_denominator);
    };
  }

  _transform_all_digital_counts(src, dst, uniform_transform);
}

/** Perform improved grey scale quantization on a color image
 *
 *  \param[in] src                 source cv::Mat of CV_8UC3
 *  \param[in] quantization_levels the number of levels to which to quantize
 *                                 the image
 *  \param[out] dst                destination cv:Mat of CV_8UC3
 */
void Igs(const cv::Mat& src, const int quantization_levels, cv::Mat& dst) {
  std::function<uint8_t(uint8_t)> igs_transform;
  float quantization_denominator = (256.0 / quantization_levels);
  // Check if it's a power of 2
  if (_is_power_of_two(quantization_levels)){
    // Determine the correct amount of bits to shift by
    uint8_t bit_shift = static_cast<uint8_t>(log2(quantization_denominator));
    uint8_t remainder = 0;
    // Create a lambda that uses bit shifting to perform IGS quantization
    // Pass the remainder by reference so it can be modified
    igs_transform = [bit_shift, &remainder](uint8_t src_digital_count) {
      // If adding the remainder would cause an overflow, make the remainder 0
      if ((remainder + src_digital_count) > 255){
        remainder = 0;
      }
      src_digital_count = src_digital_count + remainder;
      uint8_t quantized_digital_count = src_digital_count >> bit_shift;
      // Use the bitwise and operator to get the remainder
      remainder = src_digital_count & ((1 << bit_shift) - 1);
      return quantized_digital_count;
    };
  }
  else{
    float remainder = 0;
    // Create a lambda that uses division to perform IGS quantization
    // Pass the remainder by reference so it can be modified
    igs_transform = [quantization_denominator, &remainder](uint8_t src_digital_count) {
      // If adding the remainder would cause an overflow, make the remainder 0
      if ((remainder + src_digital_count) > 255){
        remainder = 0;
      }
      float src_digital_count_float = src_digital_count + remainder;
      uint8_t quantized_digital_count = static_cast<uint8_t>(src_digital_count_float / quantization_denominator);
      // Use fmod instead of modulo so it doesn't have to be an int
      remainder = fmod(src_digital_count_float, quantization_denominator);
      return quantized_digital_count;
    };
  }
  _transform_all_digital_counts(src, dst, igs_transform);
}

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
