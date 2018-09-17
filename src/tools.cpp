#include <iostream>
#include <cmath>
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
  // define return value.
  VectorXd rmse = VectorXd(4);
  rmse << 0, 0, 0, 0;
  // check input argument size.
  if ((estimations.size() != ground_truth.size()) ||
    (estimations.size() == 0)) {
    cerr << "Invalid estimations or ground_truth size." << endl;
    return rmse;
  }

  VectorXd error_square = VectorXd(4);
  error_square << 0, 0, 0, 0;
  for (unsigned int i = 0; i < estimations.size(); i++) {
    error_square = estimations[i].array() - ground_truth[i].array();
    error_square = error_square.array() * error_square.array();
    rmse += error_square;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check invalid augument input.
  if (px == 0 || py == 0) {
    cerr << "Divided by 0 occurs. This is invalid." << endl;
    return Hj;
  }

  float fact0 = px*px + py*py; // px^2 + py^2.
  float fact1 = std::sqrt(fact0); // (px^2+py^2)^1/2.
  float fact2 = std::pow(fact1, 3.0); // (px^2+py^2)^(3/2).
  // calculate jacobian matrix according to the defined formula.
  Hj << px/fact1, py/fact1, 0, 0,
    -py/fact0, px/fact0, 0, 0,
    py*(vx*py-vy*px)/fact2, px*(vy*px-vx*py)/fact2, px/fact1, py/fact1;

  return Hj;
  
}
