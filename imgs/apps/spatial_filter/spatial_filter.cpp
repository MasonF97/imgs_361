#include <ctime>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
//#include <opencv2/imgproc.hpp>

#include "imgs/ipcv/spatial_filtering/Filter2D.h"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  bool verbose = false;
  string src_filename = "";
  string dst_filename = "";

  int kernel_type = 0;

  cv::Point anchor;
  anchor.x = -1;
  anchor.y = -1;

  string border_mode_string = "constant";
  ipcv::BorderMode border_mode;
  int value = 0;

  po::options_description options("Options");
  options.add_options()("help,h", "display this message")(
      "verbose,v", po::bool_switch(&verbose), "verbose [default is silent]")(
      "source-filename,i", po::value<string>(&src_filename), "source filename")(
      "destination-filename,o", po::value<string>(&dst_filename),
      "destination filename")(
      "kernel-type,k", po::value<int>(&kernel_type),
        "kernel type (0 is blur, 1 is more blur, 2 is "
        "sharpen, 3 is Laplacian) [default is 0]")(
      "anchor-x,x", po::value<int>(&anchor.x),
      "anchor point x coord")(
      "anchor-y,y", po::value<int>(&anchor.y),
      "anchor point y coord")(
      "border-mode,m", po::value<string>(&border_mode_string),
      "border mode (constant|replicate|wrap) [default is constant]")(
      "border-value,b", po::value<int>(&value), "border value [default is 0]");

  
  po::positional_options_description positional_options;
  positional_options.add("source-filename", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
                .options(options)
                .positional(positional_options)
                .run(),
            vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << "Usage: " << argv[0] << " [options] source-filename" << endl;
    cout << options << endl;
    return EXIT_SUCCESS;
  }

  if (border_mode_string == "constant") {
    border_mode = ipcv::BorderMode::CONSTANT;
  } else if (border_mode_string == "replicate") {
    border_mode = ipcv::BorderMode::REPLICATE;
  } else if (border_mode_string == "wrap") {
    border_mode = ipcv::BorderMode::WRAP;
  } else {
    cerr << "*** ERROR *** ";
    cerr << "Provided border mode is not supported" << endl;
    return EXIT_FAILURE;
  }

  if (!boost::filesystem::exists(src_filename)) {
    cerr << "Provided source file does not exists" << endl;
    return EXIT_FAILURE;
  }

  cv::Mat src = cv::imread(src_filename, cv::IMREAD_COLOR);

  int ddepth;
  int delta;

  cv::Mat kernel;
  switch (kernel_type) {
    case 0:
      kernel.create(3, 3, CV_32FC1);
      kernel = 1;
      kernel /= 9;
      ddepth = CV_8UC3;
      delta = 0;
      break;

    case 1:
      kernel.create(5, 5, CV_32FC1);
      kernel = 1;
      kernel /= 25;
      ddepth = CV_8UC3;
      delta = 0;
      break;

    case 2:
      kernel.create(3, 3, CV_32FC1);
      kernel = -1;
      kernel.at<float>(0, 0) = 0;
      kernel.at<float>(2, 0) = 0;
      kernel.at<float>(1, 1) = 5;
      kernel.at<float>(0, 2) = 0;
      kernel.at<float>(2, 2) = 0;
      kernel /= 1;
      ddepth = CV_8UC3;
      delta = 0;
      break;

    case 3:
      kernel.create(3, 3, CV_32FC1);
      kernel = -1;
      kernel.at<float>(0, 0) = 0;
      kernel.at<float>(2, 0) = 0;
      kernel.at<float>(1, 1) = 4;
      kernel.at<float>(0, 2) = 0;
      kernel.at<float>(2, 2) = 0;
      kernel /= 1;
      ddepth = CV_8UC3;
      delta = 128;
      break;

    default:
      cerr << "*** ERROR *** ";
      cerr << "Invalid kernel type specified" << endl;
      return EXIT_FAILURE;
  }

  if (verbose) {
    cout << "Source filename: " << src_filename << endl;
    cout << "Size: " << src.size() << endl;
    cout << "Channels: " << src.channels() << endl;
    cout << "Kernel: " << endl;
    cout << kernel << endl;
    cout << "Destination filename: " << dst_filename << endl;
    cout << "Border value: " << value << endl;
    cout << "anchor_x: " << anchor.x << endl;
    cout << "anchor_y: " << anchor.y << endl;
    cout << "border mode: " << border_mode_string << endl;
  }

  cv::Mat dst;

  clock_t startTime = clock();
  uint8_t border_value = value;

  ipcv::Filter2D(src, dst, ddepth, kernel, anchor, delta, border_mode, border_value);
//  cv::filter2D(src, dst, ddepth, kernel, anchor, delta, border_type);

  clock_t endTime = clock();

  if (verbose) {
    cout << "Elapsed time: "
         << (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC)
         << " [s]" << endl;
  }

  if (dst_filename.empty()) {
    cv::imshow(src_filename, src);
    cv::imshow(src_filename + " [Filtered]", dst);
    cv::waitKey(0);
  } else {
    cv::imwrite(dst_filename, dst);
  }

  return EXIT_SUCCESS;
}
