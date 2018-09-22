// Copyright 2018 Takayuki AKAMINE
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() :
  is_initialized_(false),
  previous_timestamp_(0),
  x_state_(),
  R_laser_(),
  R_radar_(),
  H_laser_(),
  Hj_(),
  P_(),
  noise_ax(9),
  noise_ay(9)
{}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

/**
 *  Initialize vector and matrixes
 */
void FusionEKF::InitializeSelf(void) {
  // initialize state vector
  x_state_ = VectorXd(4);
  x_state_ << 1, 1, 1, 1;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
    0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;

  // process covariance matrix
  H_laser_ << 1, 0, 0, 0,
    0, 1, 0, 0;

  // Uncertainty covariance matrix
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;
}
/*
 *! Execute first step of tracking.
 @param[in] measuement_pack ; measurement dataset
 @return true: Initialization completion, false; the following tracking.
*/
bool FusionEKF::InitialTracking(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    this->InitializeSelf();
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
         Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro;
      float theta;
      ro = measurement_pack.raw_measurements_(0);
      theta = measurement_pack.raw_measurements_(1);
      x_state_(0) =  ro * cos(theta);
      x_state_(1) =  ro * sin(theta);
      ekf_.Init(x_state_, P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
         Initialize state.
      */
      x_state_(0) = measurement_pack.raw_measurements_(0);
      x_state_(1) = measurement_pack.raw_measurements_(1);
      ekf_.Init(x_state_, P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return true;
  }
  return false;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*** initialization and initial step **/
  if (InitialTracking(measurement_pack))
    // if initialization is done, return parent function.
    return;

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.setUncertainty(dt);
  ekf_.setProcessCov(dt, noise_ax, noise_ay);

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.setMeasurementNoise(R_radar_);
    ekf_.setMeasurement_Radar();
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.setMeasurementNoise(R_laser_);
    ekf_.setMeasurement_Laser(H_laser_);
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
