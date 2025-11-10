#pragma once

#include <opencv2/core.hpp>

namespace ipcv{
// Runs the fft display
void display_fourier_transform(const cv::Mat& src, bool save_video);
}


