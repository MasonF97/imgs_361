#include <ctime>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include "imgs/third_party/gnuplot-iostream.h"

#include "imgs/ipcv/otsus_threshold/OtsusThreshold.h"
#include "imgs/ipcv/utils/Utils.h"

using namespace std;

namespace po = boost::program_options;

void plot_channel_from_mat_row(const cv::Mat& pdf_mat, int row_index,
                               const std::string& title,
                               const std::string& color,
                               double threshold) {
  std::vector<std::pair<int, float>> pdf_data;
  for (int i = 0; i < pdf_mat.cols; ++i) {
      pdf_data.emplace_back(i, pdf_mat.at<float>(row_index, i));
  }

  Gnuplot gp;
  gp << "set title '" << title << " Channel PDF'\n";
  gp << "set xlabel 'Pixel Intensity'\n";
  gp << "set ylabel 'Probability Density'\n";
  gp << "set style line 1 lc rgb '" << color << "' lw 2\n";
  gp << "set style line 2 lc rgb 'black' dt 2 lw 1\n";
  gp << "plot '-' with lines ls 1 title '" << title << " PDF', \\\n";
  gp << threshold << " with lines ls 2 title 'Otsu Threshold (" << threshold << ")'\n";
  gp.send1d(pdf_data);
  gp << "pause -1\n";
}

int main(int argc, char* argv[]) {
  bool verbose = false;
  string src_filename = "";
  string dst_filename = "";

  po::options_description options("Options");
  options.add_options()("help,h", "display this message")(
      "verbose,v", po::bool_switch(&verbose), "verbose [default is silent]")(
      "source-filename,i", po::value<string>(&src_filename), "source filename")(
      "destination-filename,o", po::value<string>(&dst_filename),
      "destination filename");

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

  if (!boost::filesystem::exists(src_filename)) {
    cerr << "Provided source file does not exists" << endl;
    return EXIT_FAILURE;
  }

  cv::Mat src = cv::imread(src_filename, cv::IMREAD_COLOR);

  if (verbose) {
    cout << "Source filename: " << src_filename << endl;
    cout << "Size: " << src.size() << endl;
    cout << "Channels: " << src.channels() << endl;
    cout << "Destination filename: " << dst_filename << endl;
  }

  cv::Vec3b threshold;

  clock_t startTime = clock();

  ipcv::OtsusThreshold(src, threshold);

  clock_t endTime = clock();

  if (verbose) {
    cout << "Elapsed time: "
         << (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC)
         << " [s]" << endl;
  }

  if (verbose) {
    cout << "Threshold values = ";
    cout << threshold << endl;
  }

  if (verbose) {

    // cv::Mat_<int> dc;
    // dc = cv::Mat_<int>::zeros(3, 256);
    // for (int i = 0; i < 256; i++) {
    //   dc.at<int>(0, i) = i;
    //   dc.at<int>(1, i) = i;
    //   dc.at<int>(2, i) = i;
    // }
    cv::Mat_<int> h;
    ipcv::Histogram(src, h);
    cv::Mat_<double> pdf;
    ipcv::HistogramToPdf(h, pdf);
    // cv::Mat single_color_pdf = pdf.row(i).clone();

    plot_channel_from_mat_row(pdf, 2, "Red", "red", threshold[0]);
    plot_channel_from_mat_row(pdf, 1, "Green", "green", threshold[1]);
    plot_channel_from_mat_row(pdf, 0, "Blue", "blue", threshold[2]);

    // plot::plot2d::Params params;
    // params.set_x_label("Digital Count");
    // params.set_x_min(0);
    // params.set_x_max(255);
    // // params.xvline(threshold[0]);
    // // params.set_destination_filename(dst_filename);

    // double max_value;
    // cv::minMaxLoc(pdf, NULL, &max_value, NULL, NULL);
    // params.set_y_label("Probability Density");
    // plot::plot2d::Plot2d(dc, pdf, params);
  }

  cv::Mat lut;
  lut.create(3, 256, CV_8UC1);
  for (int b = 0; b < 3; b++) {
    for (int dc = 0; dc < 256; dc++) {
      lut.at<uint8_t>(b, dc) = (dc <= threshold[b]) ? 0 : 255;
    }
  }

  cv::Mat dst;
  ipcv::ApplyLut(src, lut, dst);

  if (dst_filename.empty()) {
    cv::imshow(src_filename, src);
    cv::imshow(src_filename + " [Thresholded]", dst);
    cv::waitKey(0);
  } else {
    cv::imwrite(dst_filename, dst);
  }

  return EXIT_SUCCESS;
}
