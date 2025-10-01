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

Eigen::Matrix3d create_rotation_matrix(double angle){
  // Create a matrix to represent a rotation transformation
  Eigen::Matrix3d rotation_matrix;
  rotation_matrix << std::cos(angle), -1 * (std::sin(angle)), 0,
                  std::sin(angle), std::cos(angle), 0,
                  0, 0, 1;
  return rotation_matrix;
}

Eigen::Matrix3d create_scale_matrix(double scale_x, double scale_y){
  // Create a matrix to represent a scale transformation
  Eigen::Matrix3d scale_matrix;
  scale_matrix << scale_x, 0, 0,
                  0, scale_y, 0,
                  0, 0, 1;
  return scale_matrix;
}

Eigen::Matrix3d create_translation_matrix(double translation_x, double translation_y){
  // Create a matrix to represent a translation transformation
  Eigen::Matrix3d translation_matrix;
  translation_matrix << 1, 0, translation_x,
                        0, 1, translation_y,
                        0, 0, 1;
  return translation_matrix;
}

Eigen::MatrixXd create_point_matrix(double x, double y){
  // Create a matrix to represent a 3x1 homogenous coordinate
  Eigen::MatrixXd dest_coords(3, 1);
  dest_coords << x,
                 y,
                 1;
  return dest_coords;
}

Eigen::Matrix3d create_transformation_matrix(double angle, double scale_x, double scale_y,
    double translation_x, double translation_y, int cols, int rows){

  // Determine how much the image needs to be translated by to move the center to the origin
  double cx = cols / 2.0;
  double cy = rows / 2.0;

  // Create a translation matrix to move the center of the image to the origin
  Eigen::Matrix3d move_center_to_origin_matrix = create_translation_matrix(-cx, -cy);
  // Create the rotate matrix
  Eigen::Matrix3d rotation_matrix = create_rotation_matrix(angle);
  // Create a translation matrix to move the center of the image back to it's original location
  Eigen::Matrix3d move_center_back_matrix = create_translation_matrix(cx, cy);
  // Create the scale matrix
  Eigen::Matrix3d scale_matrix = create_scale_matrix(scale_x, scale_y);
  // Create the translation matrix
  Eigen::Matrix3d translation_matrix = create_translation_matrix(translation_x, translation_y);

  // Combine all of the matrices into 1 affine transformation matrix
  Eigen::Matrix3d transformation_matrix = translation_matrix * scale_matrix * move_center_back_matrix * rotation_matrix * move_center_to_origin_matrix;
  return transformation_matrix;
}

vector<double> get_min_max_y(int cols, int rows, Eigen::Matrix3d transformation_matrix){
  // Define and transform each corner point
  Eigen::MatrixXd corner_1_transformed = transformation_matrix * create_point_matrix(0,0);
  Eigen::MatrixXd corner_2_transformed = transformation_matrix * create_point_matrix(cols, 0);
  Eigen::MatrixXd corner_3_transformed = transformation_matrix * create_point_matrix(0, rows);
  Eigen::MatrixXd corner_4_transformed = transformation_matrix * create_point_matrix(cols, rows);

  // Return a list with the max_x, min_x, max_y, and min_y
  vector<double> return_vector(4, 0);
  return_vector[0] =(std::max({corner_1_transformed(0,0), corner_2_transformed(0,0), corner_3_transformed(0,0), corner_4_transformed(0,0)}));
  return_vector[1] =(std::min({corner_1_transformed(0,0), corner_2_transformed(0,0), corner_3_transformed(0,0), corner_4_transformed(0,0)}));
  return_vector[2] =(std::max({corner_1_transformed(1,0), corner_2_transformed(1,0), corner_3_transformed(1,0), corner_4_transformed(1,0)}));
  return_vector[3] =(std::min({corner_1_transformed(1,0), corner_2_transformed(1,0), corner_3_transformed(1,0), corner_4_transformed(1,0)}));
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

  // Make the angle negative so that the image is rotated counter clockwise
  double ccw_angle = angle * -1;

  // Create the forward transformation matrix
  Eigen::MatrixXd transformation_matrix = create_transformation_matrix(ccw_angle, scale_x, scale_y, translation_x, translation_y, src.cols, src.rows);

  // Calculate the max/min x and y coords after the transformation to determine the output image size
  vector<double> min_max_x_y_results = get_min_max_y(src.cols, src.rows, transformation_matrix);

  double max_x = min_max_x_y_results[0];
  double min_x = min_max_x_y_results[1];
  double max_y = min_max_x_y_results[2];
  double min_y = min_max_x_y_results[3];

  // Subtract the max from the min and round it to find the output image size
  int dest_height = static_cast<int>(std::round(max_y - min_y));
  int dest_width  = static_cast<int>(std::round(max_x - min_x));

  // Create the inverse of the transformation matrix so that we can do inverse mapping
  Eigen::MatrixXd inverse_transformation_matrix = transformation_matrix.inverse();

  // Init the maps using the size information we got from transforming the corners
  map1 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  map2 = cv::Mat::zeros(dest_height, dest_width, CV_32FC1);
  // Iterate through each pixel in the destination image
  for (int y = 0; y < dest_height; y++) {
    for (int x = 0; x < dest_width; x++) {
      // Create a coordinate matrix for the pixel
      // Add the min_x and min_y to correct the location
      Eigen::MatrixXd dest_coords = create_point_matrix(x + min_x, y + min_y);
      // Multiply the inverse transformation matrix with the destination coords to determine where in the src image the pixel should come from
      Eigen::MatrixXd src_coords = inverse_transformation_matrix * dest_coords;

      map1.at<float>(y, x) = src_coords(0,0);
      map2.at<float>(y, x) = src_coords(1,0);
    }
  }
  return true;
}
}
