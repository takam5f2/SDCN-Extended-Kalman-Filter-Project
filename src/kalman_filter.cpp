#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::SetMatrixes(const float dt, const float noise_ax, const float noise_ay) {
    Tools tools;
    F_ = tools.PredictionMatrix(dt);
    Q_ = tools.CalculatePCovariance(dt, noise_ax, noise_ay);
    Hj_ = tools.CalculateJacobian(x_);
}

void KalmanFilter::Predict() {
    /**
       TODO:
       * predict the state
       */
    x_ = F_ * x_;
    MatrixXd Ftrans = F_.transpose();
    P_ = F_ * P_ * Ftrans + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
       TODO:
       * update the state by using Kalman Filter equations
       */
    VectorXd z_pred = H_ * x_;
    VectorXd error_y = z - z_pred;
    MatrixXd Htrans = H_.transpose();
    MatrixXd S = H_ * P_ * Htrans + R_; // calculate S Matrix
    MatrixXd Sinv = S.inverse();
    MatrixXd PHt = P_ * Htrans;
    MatrixXd K = PHt * Sinv; // Calculate Kalman Gain

    // update x and P
    x_ = x_ + (K * error_y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
    TODO:
    * update the state by using Extended Kalman Filter equations
    */
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    float dist_square = px*px + py*py;
    // calculate h(x')
    VectorXd H_radar = VectorXd(3);
    H_radar << std::sqrt(dist_square),
	       std::atan2(py, px),
	       (px*vx + py*vy) / std::sqrt(dist_square);
  
    VectorXd z_pred = H_radar;
    VectorXd error_y = z - z_pred;
    if (error_y(1) > M_PI) {
	error_y(1) = error_y(1) - 2*M_PI;
    }
    if (error_y(1) < -M_PI) {
	error_y(1) = error_y(1) + 2*M_PI;
    }
    
    MatrixXd Htrans = Hj_.transpose();
    MatrixXd S = Hj_ * P_ * Htrans + R_; // calculate S Matrix
    MatrixXd Sinv = S.inverse();
    MatrixXd PHt = P_ * Htrans;
    MatrixXd K = PHt * Sinv; // Calculate Kalman Gain
    // update x and P
    x_ = x_ + (K * error_y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj_) * P_;
}
