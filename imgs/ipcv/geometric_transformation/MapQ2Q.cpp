/** Implementation file for mapping a source quad on to a target quad
 *
 *  \file ipcv/geometric_transformation/MapQ2Q.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 12 Sep 2020
 */

#include "MapQ2Q.h"

#include <iostream>

#include <Eigen/Dense>
#include <opencv2/core.hpp>

using namespace std;

namespace ipcv {

/** Find the source coordinates (map1, map2) for a quad to quad mapping
 *
 *  \param[in] src       source cv::Mat of CV_8UC3
 *  \param[in] tgt       target cv::Mat of CV_8UC3
 *  \param[in] src_vertices
 *                       vertices cv:Point of the source quadrilateral (CW)
 *                       which is to be mapped to the target quadrilateral
 *  \param[in] tgt_vertices
 *                       vertices cv:Point of the target quadrilateral (CW)
 *                       into which the source quadrilateral is to be mapped
 *  \param[out] map1     cv::Mat of CV_32FC1 (size of the destination map)
 *                       containing the horizontal (x) coordinates at
 *                       which to resample the source data
 *  \param[out] map2     cv::Mat of CV_32FC1 (size of the destination map)
 *                       containing the vertical (y) coordinates at
 *                       which to resample the source data
 */
bool MapQ2Q(const cv::Mat src, const cv::Mat tgt,
            const vector<cv::Point> src_vertices,
            const vector<cv::Point> tgt_vertices, cv::Mat& map1,
            cv::Mat& map2) {

  // create the model matrix
  Eigen::MatrixXd model_matrix(8,8);
  model_matrix << tgt_vertices[0].x, tgt_vertices[0].y, 1, 0, 0, 0, tgt_vertices[0].x*-1*src_vertices[0].x, tgt_vertices[0].y*-1*src_vertices[0].x,
    tgt_vertices[1].x, tgt_vertices[1].y, 1, 0, 0, 0, tgt_vertices[1].x*-1*src_vertices[1].x, tgt_vertices[1].y*-1*src_vertices[1].x,
    tgt_vertices[2].x, tgt_vertices[2].y, 1, 0, 0, 0, tgt_vertices[2].x*-1*src_vertices[2].x, tgt_vertices[2].y*-1*src_vertices[2].x,
    tgt_vertices[3].x, tgt_vertices[3].y, 1, 0, 0, 0, tgt_vertices[3].x*-1*src_vertices[3].x, tgt_vertices[3].y*-1*src_vertices[3].x,
    0, 0, 0, tgt_vertices[0].x, tgt_vertices[0].y, 1, tgt_vertices[0].x*-1*src_vertices[0].y, tgt_vertices[0].y*-1*src_vertices[0].y,
    0, 0, 0, tgt_vertices[1].x, tgt_vertices[1].y, 1, tgt_vertices[1].x*-1*src_vertices[1].y, tgt_vertices[1].y*-1*src_vertices[1].y,
    0, 0, 0, tgt_vertices[2].x, tgt_vertices[2].y, 1, tgt_vertices[2].x*-1*src_vertices[2].y, tgt_vertices[2].y*-1*src_vertices[2].y,
    0, 0, 0, tgt_vertices[3].x, tgt_vertices[3].y, 1, tgt_vertices[3].x*-1*src_vertices[3].y, tgt_vertices[3].y*-1*src_vertices[3].y;

  // make the source coordinate vector
  Eigen::MatrixXd src_coords_matrix(8,1);
  src_coords_matrix << src_vertices[0].x,
    src_vertices[1].x,
    src_vertices[2].x,
    src_vertices[3].x,
    src_vertices[0].y,
    src_vertices[1].y,
    src_vertices[2].y,
    src_vertices[3].y;

  // Multiply the matrices to get the coefficents matrix
  Eigen::MatrixXd coefficents_matrix = model_matrix.inverse() * src_coords_matrix;
  // Use the coeffs to create the map to image matrix
  Eigen::Matrix3d map_to_image_matrix;
  map_to_image_matrix << coefficents_matrix(0,0), coefficents_matrix(1,0), coefficents_matrix(2,0),
    coefficents_matrix(3,0), coefficents_matrix(4,0), coefficents_matrix(5,0),
    coefficents_matrix(6,0), coefficents_matrix(7,0), 1;

  int dest_width = tgt.cols;
  int dest_height = tgt.rows;
  
  map1 = cv::Mat::zeros(tgt.size(), CV_32FC1);
  map2 = cv::Mat::zeros(tgt.size(), CV_32FC1);

  for(int v = 0; v < dest_height; v++){
    for (int u = 0; u < dest_width; u++){
      // Get the source coordinates
      Eigen::MatrixXd target_matrix(3,1);
      target_matrix << u,
        v,
        1;

      // Get the mapping coords by multiplying the map to image matrix and the target coords matrix
      Eigen::MatrixXd source_matrix = map_to_image_matrix * target_matrix;
      double w = source_matrix(2, 0);
      map1.at<float>(v, u) = source_matrix(0,0) / w;

      map2.at<float>(v, u) = source_matrix(1,0) / w;
      
    }
  }

  return true;
}
}
