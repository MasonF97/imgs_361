/** Interface file for a polar mapping
 *
 *  \file ipcv/geometric_transformation/MapPolar.h
 */

#pragma once

#include <iostream>

#include <Eigen/Dense>
#include <opencv2/core.hpp>

using namespace std;

namespace ipcv {

/** Find the source coordinates (map1, map2) for a polar mapping
 *
 *  \param[in] src       source cv::Mat of CV_8UC3
 *  \param[in] use_log   whether or not to do a log polar mapping
 *  \param[out] map1     cv::Mat of CV_32FC1 (size of the destination map)
 *                       containing the horizontal (x) coordinates at
 *                       which to resample the source data
 *  \param[out] map2     cv::Mat of CV_32FC1 (size of the destination map)
 *                       containing the vertical (y) coordinates at
 *                       which to resample the source data
 */
bool MapPolar(const cv::Mat src, const bool use_log, cv::Mat& map1,
            cv::Mat& map2);
}
