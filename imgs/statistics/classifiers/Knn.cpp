/** Implementation file for performing k-NN classification of a set of test
 *  images given a set of training images.
 *
 *  \file statistics/classifiers/Knn.cpp
 *  \author Carl Salvaggio, Ph.D. (salvaggio@cis.rit.edu)
 *  \date 22 Nov 2023
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

using namespace std;
#include <opencv2/core.hpp>
#include <algorithm>
#include <utility>

#include "imgs/statistics/classifiers/Knn.h"
#include "imgs/statistics/data_readers/Mnist.h"

namespace statistics {


double minkowski_distance(double p, std::vector<float>& vector_1, std::vector<float>& vector_2){
  size_t n = vector_1.size();
  double sum = 0.0;
  // Check if p is 1 or 2, if so use a more optimized calculation
  if(p == 1){
    return cv::norm(vector_1, vector_2, cv::NORM_L1);
  }else if(p == 2){
    return cv::norm(vector_1, vector_2, cv::NORM_L2);
  }else{
    // Regular minkowski formula
    for (size_t i = 0; i < n; i++){
      sum += std::pow(std::abs(vector_1[i] - vector_2[i]), p);
    }
    return std::pow(sum, 1.0 / p);
  }
}

std::vector<int> find_k_closest(int k, std::vector<double>& distances){
  // Create an indexes array with all of the indexes
  std::vector<int> indexes(distances.size());
  std::iota(indexes.begin(), indexes.end(), 0);
  // Sort the indexes array based on the distances array
  std::nth_element(indexes.begin(), indexes.begin() + k, indexes.end(),
    [&](int a, int b){ return distances[a] < distances[b]; });
  // Only return k elements
  indexes.resize(k);
  return indexes;
}

unsigned char predict_label(const std::vector<unsigned char>& training_labels, std::vector<int>& k_closest){
  // Create a list that stores the amount of times a label is present
  int label_counts[10] = {0};
  for (int x : k_closest){
    label_counts[training_labels[x]]++;
  }
  // Find and return the label that has the highest count
  int highest_count_label = 0;
  int highest_count = 0;
  for (int i = 0; i < 10; i++) {
    if (label_counts[i] > highest_count) {
      highest_count = label_counts[i];
      highest_count_label = i;
    }
  }
  return static_cast<unsigned char>(highest_count_label);
}

std::vector<unsigned char> Knn(
    const std::vector<cv::Mat>& test_images,
    const std::vector<cv::Mat>& training_images,
    const std::vector<unsigned char>& training_labels, const int k,
    const double p) {
  // Instantiate a vector to hold the predicted label for each test image
  std::vector<unsigned char> predicted_test_labels;

  int n_train_size = training_images.size();
  int n_test_size = test_images.size();

  // transform all of the training images to vectors
  std::vector<std::vector<float>> training_vectors;
  training_vectors.reserve(training_images.size());
  for (const cv::Mat& train_img : training_images) {
    // convert to float
    cv::Mat img_float;
    train_img.convertTo(img_float, CV_32F, 1.0 / 255.0);

    std::vector<float> vec(img_float.total());
    std::memcpy(vec.data(), img_float.ptr<float>(), img_float.total() * sizeof(float));

    training_vectors.push_back(std::move(vec));
  }

  // transform all of the test images to vectors
  std::vector<std::vector<float>> test_vectors;
  test_vectors.reserve(test_images.size());
  for (const cv::Mat& test_img : test_images) {
    // convert to float
    cv::Mat img_float;
    test_img.convertTo(img_float, CV_32F, 1.0 / 255.0);

    std::vector<float> vec(img_float.total());
    std::memcpy(vec.data(), img_float.ptr<float>(), img_float.total() * sizeof(float));

    test_vectors.push_back(std::move(vec));
  }

  // KNN
  // label each test image
  std::vector<double> distances(training_images.size());
  for(int i = 0;i<n_test_size;i++){
    // Find the distances to all of the training vectors
    for(int t = 0;t<n_train_size;t++){
      distances[t] = minkowski_distance(p, test_vectors[i], training_vectors[t]);
    }
    // find the k closest points to the test point
    std::vector<int> k_closest = find_k_closest(k, distances);
    // Predict a label based on the k_closest
    predicted_test_labels.push_back(predict_label(training_labels, k_closest));
  }
  return predicted_test_labels;
}
}
