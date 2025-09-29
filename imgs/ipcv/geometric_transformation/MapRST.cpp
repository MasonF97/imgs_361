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

// temp
// template <typename T>
// void printMatrix(const std::vector<std::vector<T>>& matrix) {
//   for (const auto& row : matrix) {
//     for (const auto& val : row) {
//       std::cout << val << " ";
//     }
//     std::cout << std::endl;
//   }
// }

vector<vector<double>> multiply_matrices(std::vector<std::vector<double>> left_matrix, vector<vector<double>> right_matrix){
  int left_rows = static_cast<int>(left_matrix.size());
  int left_cols = static_cast<int>(left_matrix[0].size());
  int right_rows = static_cast<int>(right_matrix.size());
  int right_cols = static_cast<int>(right_matrix[0].size());

  if (left_cols != right_rows){
    throw std::invalid_argument("Matrices are not able to be multiplied.");
  }
  vector<vector<double>> output_matrix(left_rows, std::vector<double>(right_cols, 0));
  for (int i = 0; i <left_rows; i++ ){
    for (int j = 0; j < right_cols; j++){
      for (int k = 0; k < left_cols; k++){
        output_matrix[i][j] += left_matrix[i][k] * right_matrix[k][j];
      }
    }
  }
  return output_matrix;
}

vector<vector<double>> create_rotation_matrix(double angle, bool inverse = false){
  vector<vector<double>> scale_matrix(3, vector<double>(3, 0));
  if (inverse){
    scale_matrix[0][1] = std::sin(angle);
    scale_matrix[1][0] = -1 * (std::sin(angle));
  }else{
    scale_matrix[0][1] = -1 * (std::sin(angle));
    scale_matrix[1][0] = std::sin(angle);
  }
  scale_matrix[0][0] = std::cos(angle);
  scale_matrix[1][1] = std::cos(angle);
  scale_matrix[2][2] = 1;
  return scale_matrix;
}

vector<vector<double>> create_scale_matrix(double scale_x, double scale_y, bool inverse = false){
  vector<vector<double>> scale_matrix(3, vector<double>(3, 0));
  if (inverse){
    scale_matrix[0][0] = 1 / scale_x;
    scale_matrix[1][1] = 1 / scale_y;
  }else{
    scale_matrix[0][0] = scale_x;
    scale_matrix[1][1] = scale_y;
  }
  scale_matrix[2][2] = 1;
  return scale_matrix;
}

vector<vector<double>> create_translation_matrix(double translation_x, double translation_y, bool inverse = false){
  vector<vector<double>> scale_matrix(3, vector<double>(3, 0));
  if (inverse){
    scale_matrix[0][2] = -1 * translation_x;
    scale_matrix[1][2] = -1 * translation_y;
  }else{
    scale_matrix[0][2] = translation_x;
    scale_matrix[1][2] = translation_y;
  }
  scale_matrix[0][0] = 1;
  scale_matrix[1][1] = 1;
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

vector<vector<double>> create_transformation_matrix(double angle, double scale_x, double scale_y,
    double translation_x, double translation_y, int cols, int rows, bool inverse = false){

  double cx = cols / 2.0;
  double cy = rows / 2.0;

  vector<vector<double>> move_center_to_origin_matrix = create_translation_matrix(-cx, -cy, inverse);
  
  // Create the rotate matrix
  vector<vector<double>> rotation_matrix = create_rotation_matrix(angle, inverse);

  vector<vector<double>> move_center_back_matrix = create_translation_matrix(cx, cy, inverse);

  // Create the scale matrix
  vector<vector<double>> scale_matrix = create_scale_matrix(scale_x, scale_y, inverse);

  // Create the translation matrix
  vector<vector<double>> translation_matrix = create_translation_matrix(translation_x, translation_y, inverse);

  if (inverse){
    return multiply_matrices(move_center_to_origin_matrix, 
      multiply_matrices(rotation_matrix,
        multiply_matrices(move_center_back_matrix,
          multiply_matrices(scale_matrix, translation_matrix))));
  }else{
    return multiply_matrices(translation_matrix, 
      multiply_matrices(scale_matrix,
        multiply_matrices(move_center_back_matrix,
          multiply_matrices(rotation_matrix, move_center_to_origin_matrix))));
  }
}

vector<double> get_min_max_y(int cols, int rows, vector<vector<double>> transformation_matrix){
  vector<vector<double>> corner_1 = create_point_matrix(0,0);
  vector<vector<double>> corner_2 = create_point_matrix(cols, 0);
  vector<vector<double>> corner_3 = create_point_matrix(0, rows);
  vector<vector<double>> corner_4 = create_point_matrix(cols, rows);

  vector<vector<double>> corner_1_transformed = multiply_matrices(transformation_matrix, corner_1);
  vector<vector<double>> corner_2_transformed = multiply_matrices(transformation_matrix, corner_2);
  vector<vector<double>> corner_3_transformed = multiply_matrices(transformation_matrix, corner_3);
  vector<vector<double>> corner_4_transformed = multiply_matrices(transformation_matrix, corner_4);

  cout << corner_1_transformed[0][0] << " " << corner_1_transformed[1][0] << endl;
  cout << corner_2_transformed[0][0] << " " << corner_2_transformed[1][0] << endl;
  cout << corner_3_transformed[0][0] << " " << corner_3_transformed[1][0] << endl;
  cout << corner_4_transformed[0][0] << " " << corner_4_transformed[1][0] << endl;

  vector<double> return_vector(4, 0);
  return_vector[0] =(std::max({corner_1_transformed[0][0], corner_2_transformed[0][0], corner_3_transformed[0][0], corner_4_transformed[0][0]}));
  return_vector[1] =(std::min({corner_1_transformed[0][0], corner_2_transformed[0][0], corner_3_transformed[0][0], corner_4_transformed[0][0]}));
  return_vector[2] =(std::max({corner_1_transformed[1][0], corner_2_transformed[1][0], corner_3_transformed[1][0], corner_4_transformed[1][0]}));
  return_vector[3] =(std::min({corner_1_transformed[1][0], corner_2_transformed[1][0], corner_3_transformed[1][0], corner_4_transformed[1][0]}));
  return return_vector;
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

  double ccw_angle = angle * -1;

  vector<vector<double>> transformation_matrix = create_transformation_matrix(ccw_angle, scale_x, scale_y, translation_x, translation_y, src.cols, src.rows);

  vector<double> min_max_x_y_results = get_min_max_y(src.cols, src.rows, transformation_matrix);

  double max_x = min_max_x_y_results[0];
  double min_x = min_max_x_y_results[1];
  double max_y = min_max_x_y_results[2];
  double min_y = min_max_x_y_results[3];

  cout << endl;
  cout << max_x << endl;
  cout << min_x << endl;
  cout << max_y << endl;
  cout << min_y << endl;

  int dest_height = max_y - min_y;
  int dest_width = max_x - min_x;

  cout << dest_height << endl;
  cout << dest_width << endl;

  vector<vector<double>> inverse_transformation_matrix = create_transformation_matrix(ccw_angle, scale_x, scale_y, translation_x, translation_y,src.cols, src.rows, true);

  map1 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  map2 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  for (int y = 0; y < dest_height; y++) {
    for (int x = 0; x < dest_width; x++) {
      vector<vector<double>> dest_coords(3, std::vector<double>(1, 0));
      dest_coords[0][0] = x + min_x;
      dest_coords[1][0] = y + min_y;
      dest_coords[2][0] = 1;
      vector<vector<double>> src_coords = multiply_matrices(inverse_transformation_matrix, dest_coords);
      
      // if (src_coords[0][0] < 0){
      //   cout << src_coords[0][0] << endl;
      // }
      // if (src_coords[1][0] < 0){
      //   cout << src_coords[1][0] << endl;
      // }

      map1.at<float>(y, x) = src_coords[0][0];
      map2.at<float>(y, x) = src_coords[1][0];
    }
  }
  return true;
}
}
