/** Implementation file for finding map coordinates for an RST transformation
 *
 *  \file ipcv/geometric_transformation/MapRST.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 26 Sep 2019
 */

#include "MapRST.h"

#include <iostream>

#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

using namespace std;

namespace ipcv {

vector<vector<double>> multiply_matrices(std::vector<std::vector<double>> left_matrix, vector<vector<double>> right_matrix){
  int left_rows = static_cast<int>(left_matrix.size());
  int left_cols = static_cast<int>(left_matrix[0].size());
  int right_rows = static_cast<int>(right_matrix.size());
  int right_cols = static_cast<int>(right_matrix[0].size());

  if (left_cols == right_rows){
    throw std::invalid_argument("Matrices are not able to be multiplied.");
  }
  vector<vector<double>> output_matrix(left_rows, std::vector<double>(right_cols, 0));
  for (int i = 0; i <left_rows; i++ ){
    for (int j = 0; i < right_cols; i++){
      for (int k = 0; k < left_cols; k++){
        output_matrix[i][j] += left_matrix[i][k] * right_matrix[k][j];
      }
    }
  }
  return output_matrix;
}

vector<vector<double>> create_rotation_matrix(double angle){
  vector<vector<double>> scale_matrix(3, vector<double>(3, 0));
  scale_matrix[0][0] = std::cos(angle);
  scale_matrix[0][1] = -1 * (std::sin(angle));
  scale_matrix[1][0] = std::sin(angle);
  scale_matrix[1][1] = std::cos(angle);
  scale_matrix[2][2] = 1;
  return scale_matrix;
}

vector<vector<double>> create_scale_matrix(double scale_x, double scale_y){
  vector<vector<double>> scale_matrix(3, vector<double>(3, 0));
  scale_matrix[0][0] = scale_x;
  scale_matrix[1][1] = scale_y;
  scale_matrix[2][2] = 1;
  return scale_matrix;
}

vector<vector<double>> create_translation_matrix(double translation_x, double translation_y){
  vector<vector<double>> scale_matrix(3, vector<double>(3, 0));
  scale_matrix[0][0] = 1;
  scale_matrix[1][1] = 1;
  scale_matrix[0][2] = translation_x;
  scale_matrix[1][2] = translation_y;
  scale_matrix[2][2] = 1;
  return scale_matrix;
}

vector<vector<double>> create_point_matrix(double x, double y){
  vector<vector<double>> dest_coords(3, std::vector<double>(1, 0));
  dest_coords[0][0] = x;
  dest_coords[1][0] = y;
  dest_coords[2][0] = 1;
  return dest_coords;
}


/** Find the map coordinates (map1, map2) for an RST transformation
 *
 *  \param[in] src           source cv::Mat of CV_8UC3
 *  \param[in] angle         rotation angle (CCW) [radians]
 *  \param[in] scale_x       horizontal scale
 *  \param[in] scale_y       vertical scale
 *  \param[in] translation_x horizontal translation [+ right]
 *  \param[in] translation_y vertical translation [+ up
 *  \param[out] map1         cv::Mat of CV_32FC1 (size of the destination map)
 *                           containing the horizontal (x) coordinates at
 *                           which to resample the source data
 *  \param[out] map2         cv::Mat of CV_32FC1 (size of the destination map)
 *                           containing the vertical (y) coordinates at
 *                           which to resample the source data
 */
bool MapRST(const cv::Mat src, const double angle, const double scale_x,
            const double scale_y, const double translation_x,
            const double translation_y, cv::Mat& map1, cv::Mat& map2) {

  // Create transformation matrices and then combine them

  // Create the rotate matrix
  vector<vector<double>> rotation_matrix = create_rotation_matrix(angle);

  // Create the scale matrix
  vector<vector<double>> scale_matrix = create_scale_matrix(scale_x, scale_y);

  // Create the translation matrix
  vector<vector<double>> translation_matrix = create_translation_matrix(scale_x, scale_y);


  vector<vector<double>> scale_and_translate_matrix = multiply_matrices(translation_matrix, scale_matrix);
  vector<vector<double>> transformation_matrix = multiply_matrices(scale_and_translate_matrix, rotation_matrix);

  vector<vector<double>> corner1 = create_point_matrix(0,0);
  vector<vector<double>> corner2 = create_point_matrix(0,src.cols);
  vector<vector<double>> corner3 = create_point_matrix(src.rows,0);
  vector<vector<double>> corner4 = create_point_matrix(src.rows,src.cols);


  vector<vector<double>> corner_1_transformed = multiply_matrices(transformation_matrix, corner1);
  vector<vector<double>> corner_2_transformed = multiply_matrices(transformation_matrix, corner2);
  vector<vector<double>> corner_3_transformed = multiply_matrices(transformation_matrix, corner3);
  vector<vector<double>> corner_4_transformed = multiply_matrices(transformation_matrix, corner4);

  double max_x = std::min({corner_1_transformed[0][0], corner_2_transformed[0][0], corner_3_transformed[0][0], corner_4_transformed[0][0]});
  double min_x = std::min({corner_1_transformed[0][0], corner_2_transformed[0][0], corner_3_transformed[0][0], corner_4_transformed[0][0]});
  double max_y = std::min({corner_1_transformed[0][1], corner_2_transformed[0][1], corner_3_transformed[0][1], corner_4_transformed[0][1]});
  double min_y = std::min({corner_1_transformed[0][1], corner_2_transformed[0][1], corner_3_transformed[0][1], corner_4_transformed[0][1]});

  int dest_width = max_x - min_x;
  int dest_height = max_y - min_y;

  cout << dest_height << endl;
  cout << dest_width << endl;

  // map1 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  // map2 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  // for (int y = 0; y < dest_height; y++) {
  //   for (int x = 0; x < dest_width; x++) {
  //     vector<vector<double>> dest_coords(3, std::vector<double>(1, 0));
  //     dest_coords[0][0] = x;
  //     dest_coords[1][0] = y;
  //     dest_coords[2][0] = 1;
  //     vector<vector<double>> src_coords = multiply_matrices(transformation_matrix, dest_coords);
      
  //     map1.at<float>(y, x) = src_coords[0][0];
  //     map2.at<float>(y, x) = src_coords[1][0];
  //   }
  // }



  return true;
}
}
