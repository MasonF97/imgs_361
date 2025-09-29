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

cv::Vec3f bilinear_interpolation(vector<cv::Vec3b> point_dcs, float x, float y, int dc_i, int dc_j){
  float dx = x - dc_j;
  float dy = y - dc_i;

  
  cv::Vec3f dc_x_j = (point_dcs[1] - point_dcs[0])* dx + point_dcs[0];
  cv::Vec3f dc_x_j1 = (point_dcs[3] - point_dcs[2])* dx + point_dcs[2];

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

  for (int i = 0; i < dst.rows;i++){
    for (int j = 0; j < dst.cols; j++){
      float x = map1.at<float>(i, j);
      float y = map2.at<float>(i, j);
      if (interpolation == Interpolation::NEAREST){
        int src_x = static_cast<int>(x + 0.5);
        int src_y = static_cast<int>(y + 0.5);
        if (src_x < 0 || src_y < 0 || src_x >= src.cols || src_y >= src.rows){
          if (border_mode == BorderMode::CONSTANT){
            dst.at<cv::Vec3b>(i, j) = cv::Vec3b (border_value, border_value, border_value);
          }else{
            int clamped_y = std::clamp(src_y, 0, src.rows-1);
            int clamped_x = std::clamp(src_x, 0, src.cols-1);
            dst.at<cv::Vec3b>(i, j) = src.at<cv::Vec3b>(clamped_y, clamped_x);
          }
        }else{
          dst.at<cv::Vec3b>(i, j) = src.at<cv::Vec3b>(src_y, src_x);
        }
      }else{
        int dc_i = std::floor(y);
        int dc_j = std::floor(x);

        vector<vector<int>> dc_point_coords = {
            {dc_i, dc_j},
            {dc_i + 1, dc_j},
            {dc_i, dc_j + 1},
            {dc_i + 1, dc_j + 1}
          };

        vector<cv::Vec3b> dc_points_dc(4);

        for( int k = 0; k< 4; k ++){
          vector<int> point = dc_point_coords[k];
          if (point[0] < 0 || point[1] < 0 || point[1] >= src.cols || point[0] >= src.rows){
            if (border_mode == BorderMode::CONSTANT){
              dc_points_dc[k] = cv::Vec3b (border_value, border_value, border_value);
            }else{
              int clamped_y = clamp(point[0], 0, src.rows-1);
              int clamped_x = clamp(point[1], 0, src.cols-1);
              dc_points_dc[k] = src.at<cv::Vec3b>(clamped_y, clamped_x);
            }
          }else{
            dc_points_dc[k] = src.at<cv::Vec3b>(point[0], point[1]);
          }
        }

        dst.at<cv::Vec3b>(i, j) = bilinear_interpolation(dc_points_dc, x, y, dc_i, dc_j);
      }
    }
  }

  return true;
}
}
