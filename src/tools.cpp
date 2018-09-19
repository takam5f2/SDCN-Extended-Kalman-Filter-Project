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
    MatrixXd Hj_ = MatrixXd(3, 4);
    
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    float c1 = px*px + py*py;
    float c2 = sqrt(c1);
    float c3 = c1 * c2;
    if(fabs(c1) < 0.0001){
	cout << "CalculateJacobian () - Error - Division by Zero" << endl;
	return Hj_;
    } 
    Hj_ << (px/c2), (py/c2), 0, 0,
	-(py/c1), (px/c1), 0, 0,
	py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    return Hj_;
  
}

MatrixXd Tools::PredictionMatrix(const float delta_t) {
    MatrixXd F_;
    F_ = MatrixXd(4,4);
    F_ << 1, 0, delta_t, 0,
	       0, 1, 0, delta_t,
	       0, 0, 1, 0,
	       0, 0, 0, 1;
    return F_;
}

MatrixXd Tools::CalculatePCovariance(const float delta_t, const float noise_ax, const float noise_ay) {
    MatrixXd Q_ ;
    float dt_2 = delta_t * delta_t;
    float dt_3 = dt_2 * delta_t;
    float dt_4 = dt_3 * delta_t;
    Q_ = MatrixXd(4,4);
    Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
               0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
	       dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
	       0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
    return Q_;
}
