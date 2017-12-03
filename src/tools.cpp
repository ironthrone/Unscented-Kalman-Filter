#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  VectorXd se(4);
  se << 0, 0, 0, 0;
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd different = (estimations[i] - ground_truth[i]);
    different = different.array() * different.array();
    se += different;
  }

  rmse = (se / estimations.size());
  rmse = rmse.array().sqrt();

  return rmse;
}