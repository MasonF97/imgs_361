#include "fourier_transform_display.h"

#include <iostream>
#include "imgs/ipcv/utils/Utils.h"

using namespace std;

namespace ipcv {

void display_fourier_transform(const cv::Mat& src, bool save_video) {
    // prepare and compute the FFT using the pipeline shown in the example in class
    cv::Mat spectrum = src | pcv::to_gray() | pcv::to_double() | pcv::pad_to_even() | pcv::dft();
    cv::Mat log_mag = spectrum | pcv::magnitude() | pcv::log() | pcv::fft_shift() | pcv::normalize();

    // Create the mats that will be filled as we add components
    cv::Mat used_coeffs = cv::Mat::zeros(spectrum.size(), CV_64FC2);
    cv::Mat summed = cv::Mat::zeros(src.size(), CV_64F);
    cv::Mat current_component;
    cv::Mat current_component_scaled;

    bool paused = false;

    // create a named window for the display
    cv::namedWindow("(Right to Left, Top to Bottom) Original | Fourier Transform log(mag) | Coeffs Used log(mag) | Current Component | Current Component Scaled | Summed Components");

    // Setup the video writer
    std::string output_filename = "../imgs/ipcv/fourier_transform/fft_display.avi";
    // 'm','p','4','v' wouldn't work, so I had to use m j p g
    int codec = cv::VideoWriter::fourcc('M','J','P','G');
    double fps = 30.0;
    cv::Size frame_size(src.cols*3, src.rows*2);

    // Open the video writer
    cv::VideoWriter writer(output_filename, codec, fps, frame_size,false);
    if (!writer.isOpened()) {
        std::cerr << "Error: Could not open the output video for write\n";
        return;
    }

    // Combine src and the log mag for the display
    // Need to convert them to greyscale and 64f because everything needs to be the same type
    cv::Mat src_and_log_mag;
    cv::Mat src_64F;
    (src | pcv::to_gray()).convertTo(src_64F, CV_64F);
    cv::hconcat(src_64F | pcv::normalize(), log_mag | pcv::normalize(), src_and_log_mag);

    // setup ints for the loop
    int total_rows = spectrum.rows;
    int total_cols = spectrum.cols;
    int step = 0;

    cout << "Press 'p' to pause/resume, 'q' or ESC to quit.\n";

    while (true) {
        if (!paused) {
            // get the next frequency coordinate to add
            int r = step / total_cols;
            int c = step % total_cols;
            if (r >= total_rows){
                r = total_rows - 1;
            }

            // create the zero matrix except current frequency component
            cv::Mat single_freq = cv::Mat::zeros(spectrum.size(), CV_64FC2);
            single_freq.at<cv::Vec2d>(r, c) = spectrum.at<cv::Vec2d>(r, c);

            // use inverse dft to get the current component
            current_component = single_freq |
                pcv::dft(cv::DFT_INVERSE + cv::DFT_SCALE + cv::DFT_REAL_OUTPUT);

            // normalize the current component so we can actually see it.
            current_component_scaled = current_component | pcv::normalize();

            // add the current component to the running summed image
            summed += current_component;

            // mark coeff as used
            used_coeffs.at<cv::Vec2d>(r, c) = spectrum.at<cv::Vec2d>(r, c);

            // I concated all of the displays as suggested by the prof

            // hconcat requires the name number of rows and vconcat requires the same amount of columns
            // Combine the first 3 horizontally
            cv::Mat first_3;
            cv::hconcat(src_and_log_mag, used_coeffs | pcv::magnitude() | pcv::log() | pcv::fft_shift() | pcv::normalize(), first_3);            

            // Combine the last 3 horizontally
            cv::Mat fourth_and_fifth;
            cv::hconcat(current_component, current_component_scaled, fourth_and_fifth);
            cv::Mat last_3;
            cv::hconcat(fourth_and_fifth, summed | pcv::normalize(), last_3);

            // Merge them vertically
            cv::Mat merged_display;
            cv::vconcat(first_3, last_3, merged_display);

            // Show the merged displays
            cv::imshow("(Right to Left, Top to Bottom) Original | Fourier Transform log(mag) | Coeffs Used log(mag) | Current Component | Current Component Scaled | Summed Components", merged_display);

            // If we are saving the video, write the display
            if(save_video){
                // We have to convert it to cv_8U to save it
                cv::Mat merged_display_8U;
                merged_display.convertTo(merged_display_8U, CV_8U, 255.0);
                writer.write(merged_display_8U);
            }

            step++;
            // check if all of the components have been added
            if (step >= total_rows * total_cols){
                // if they have been, stop the loop
                break;
            }
        }

        // if they press p, pause it
        // if they press 1, quit
        int key = cv::waitKey(10);
        if (key == 'p' || key == 'P'){
            paused = !paused;
        }
        else if (key == 27 || key == 'q' || key == 'Q'){
            break;
        }
    }
    writer.release();
    cv::destroyAllWindows();
}

}