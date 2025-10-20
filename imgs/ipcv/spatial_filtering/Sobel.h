/** Interface file for a sobel gradient
 *
 *  \file ipcv/spatial_filtering/Sobel.h
 */

#pragma once

#include <opencv2/core.hpp>
#include "Filter2D.h"

namespace ipcv {

/** Does a Sobel gradient
 *
 *  \param[in] src          source cv::Mat of CV_8UC3
 *  \param[out] dst         destination cv::Mat of ddepth type
 *  \param[in] delta        optional value added to the filtered pixels
 *                          before storing them in dst
 *  \param[in] border_mode  pixel extrapolation method
 *  \param[in] border_value value to use for constant border mode
 */
bool Sobel(const cv::Mat& src, cv::Mat& dst,
                 const int delta = 0,
                 const BorderMode border_mode = BorderMode::WRAP,
                 uint8_t border_value = 0);
}
