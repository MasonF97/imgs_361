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
    for (size_t i = 0; i < n; i++) {
      sum += std::abs(vector_1[i] - vector_2[i]);
    }
    return sum;
  }else if(p == 2){
    for (size_t i = 0; i < n; i++) {
      double diff = vector_1[i] - vector_2[i];
      sum += diff * diff;
    }
    return std::sqrt(sum);
  }else{
    // Regular minkowski formula
    for (size_t i = 0; i < n; i++){
      sum += std::pow(std::abs(vector_1[i] - vector_2[i]), p);
    }
    return std::pow(sum, 1.0 / p);
  }
}

std::vector<int> find_k_closest(int k, std::vector<double>& distances){
  // make an array of pairs so we know the indexes
  std::vector<std::pair<double, int>> distances_with_indexes;
  distances_with_indexes.reserve(distances.size());
  for(size_t i = 0; i < distances.size(); i++){
    distances_with_indexes.push_back({distances[i], (int)i});
  }
  // sort the distances
  std::sort(distances_with_indexes.begin(), distances_with_indexes.end());

  // get the first k values
  std::vector<int> k_closest;
  for(int i = 0; i<k;i++){
    k_closest.push_back(distances_with_indexes[i].second);
  }
  return k_closest;
}

unsigned char predict_label(const std::vector<unsigned char>& training_labels, std::vector<int>& k_closest){
  // Get the labels of the k closest
  std::vector<unsigned char> k_closest_labels;
  for (size_t i = 0;i<k_closest.size();i++){
    k_closest_labels.push_back(training_labels[k_closest[(int)i]]);
  }
  // Create a map that stores the amount of times a label is present
  std::map<unsigned char, int> label_counts;
  for (unsigned char x : k_closest_labels){
    label_counts[x]++;
  }
  // Find and return the label that has the highest count
  unsigned char highest_count_label = k_closest_labels[0];
  int highest_count = 0;
  for (const auto& pair : label_counts) {
    if (pair.second > highest_count) {
      highest_count = pair.second;
      highest_count_label = pair.first;
    }
  }
  return highest_count_label;
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

  int counter = 0;

  // KNN
  // label each test image
  std::vector<double> distances;
  for(int i = 0;i<n_test_size;i++){
    // Find the distances to all of the training vectors
    for(int t = 0;t<n_train_size;t++){
      distances.push_back(minkowski_distance(p, test_vectors[i], training_vectors[t]));
    }
    // find the k closest points to the test point
    std::vector<int> k_closest = find_k_closest(k, distances);
    // Predict a label based on the k_closest
    predicted_test_labels.push_back(predict_label(training_labels, k_closest));
    counter++;
    if(counter % 100 == 0){
      cout << counter << endl;
    }
    // cout << counter << endl;
  }
  return predicted_test_labels;
}
}
