/** Implementation file for remapping source values to map locations
 *
 *  \file ipcv/geometric_transformation/Remap.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 15 Sep 2018
 */

#include "Remap.h"

#include <iostream>
#include <algorithm>

using namespace std;

namespace ipcv {

cv::Vec3f bilinear_interpolation(vector<cv::Vec3b> neighbor_pixels_dc, float x, float y, int neighbor_i, int neighbor_j){
  // Perform bilinear interpolation for the x,y coordinate using four surrounding pixels
  float dx = x - neighbor_j;
  float dy = y - neighbor_i;

  // Interpolate along x for top and bottom rows
  cv::Vec3f dc_x_j = (neighbor_pixels_dc[1] - neighbor_pixels_dc[0])* dx + neighbor_pixels_dc[0];
  cv::Vec3f dc_x_j1 = (neighbor_pixels_dc[3] - neighbor_pixels_dc[2])* dx + neighbor_pixels_dc[2];

  // Interpolate along y between the two rows
  return (dc_x_j1 - dc_x_j) * dy + dc_x_j;
}

/** Remap source values to the destination array at map1, map2 locations
 *
 *  \param[in] src            source cv::Mat of CV_8UC3
 *  \param[out] dst           destination cv::Mat of CV_8UC3 for remapped values
 *  \param[in] map1           cv::Mat of CV_32FC1 (size of the destination map)
 *                            containing the horizontal (x) coordinates at
 *                            which to resample the source data
 *  \param[in] map2           cv::Mat of CV_32FC1 (size of the destination map)
 *                            containing the vertical (y) coordinates at
 *                            which to resample the source data
 *  \param[in] interpolation  interpolation to be used for resampling
 *  \param[in] border_mode    border mode to be used for out of bounds pixels
 *  \param[in] border_value   border value to be used when constant border mode
 *                            is to be used
 */
bool Remap(const cv::Mat& src, cv::Mat& dst, const cv::Mat& map1,
           const cv::Mat& map2, const Interpolation interpolation,
           const BorderMode border_mode, const uint8_t border_value) {
  dst.create(map1.size(), src.type());

  // Loop over each pixel in the destination image
  for (int i = 0; i < dst.rows;i++){
    for (int j = 0; j < dst.cols; j++){
      // Get the remap coordinates from the maps
      float x = map1.at<float>(i, j);
      float y = map2.at<float>(i, j);
      if (interpolation == Interpolation::NEAREST){
        // Nearest neighbor interpolation
        int src_x = static_cast<int>(x + 0.5);
        int src_y = static_cast<int>(y + 0.5);
        // Check if any of the pixels are out of bounds
        if (src_x < 0 || src_y < 0 || src_x >= src.cols || src_y >= src.rows){
          if (border_mode == BorderMode::CONSTANT){
            // If border mode is constant, set the pixel to the border_value
            dst.at<cv::Vec3b>(i, j) = cv::Vec3b (border_value, border_value, border_value);
          }else{
            // If border mdoe is replicate, clamp the coords to get the nearest pixel in the image
            int clamped_y = std::clamp(src_y, 0, src.rows-1);
            int clamped_x = std::clamp(src_x, 0, src.cols-1);
            dst.at<cv::Vec3b>(i, j) = src.at<cv::Vec3b>(clamped_y, clamped_x);
          }
        }else{
          dst.at<cv::Vec3b>(i, j) = src.at<cv::Vec3b>(src_y, src_x);
        }
      }else{
        // find two of the neighbors for bilinear interpolation
        int neighbor_i = std::floor(y);
        int neighbor_j = std::floor(x);

        // Create a vector to store the neighbor pixel coords
        vector<vector<int>> neighbor_pixel_coords = {
            {neighbor_i, neighbor_j},
            {neighbor_i + 1, neighbor_j},
            {neighbor_i, neighbor_j + 1},
            {neighbor_i + 1, neighbor_j + 1}
          };

        // Create a vector to store the neighbor pixels digital counts
        vector<cv::Vec3b> neighbor_pixels_dc(4);

        for( int k = 0; k< 4; k ++){
          // For each of the neighbors, check if they are out of bounds
          vector<int> point = neighbor_pixel_coords[k];
          if (point[0] < 0 || point[1] < 0 || point[1] >= src.cols || point[0] >= src.rows){
            if (border_mode == BorderMode::CONSTANT){
              // If the border mode is constant, set that pixels value to the border value
              neighbor_pixels_dc[k] = cv::Vec3b (border_value, border_value, border_value);
            }else{
              // If the border mode is replicate, clamp the coords to find the closest pixel in the image
              int clamped_y = clamp(point[0], 0, src.rows-1);
              int clamped_x = clamp(point[1], 0, src.cols-1);
              neighbor_pixels_dc[k] = src.at<cv::Vec3b>(clamped_y, clamped_x);
            }
          }else{
            neighbor_pixels_dc[k] = src.at<cv::Vec3b>(point[0], point[1]);
          }
        }
        // Use bilinear interpolation to find the new value for the destination pixel
        dst.at<cv::Vec3b>(i, j) = bilinear_interpolation(neighbor_pixels_dc, x, y, neighbor_i, neighbor_j);
      }
    }
  }

  return true;
}
}
