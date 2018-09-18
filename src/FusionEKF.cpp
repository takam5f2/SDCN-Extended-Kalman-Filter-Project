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
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
	
    H_laser_ << 1, 0, 0, 0,
	        0, 1, 0, 0;

    /**
       TODO:
       * Finish initializing the FusionEKF.
       * Set the process and measurement noises
       */
    P_ = MatrixXd(4,4);
    P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;
    noise_ax = 9;
    noise_ay = 9;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    // Eigen::MatrixXd *R_;
    // Eigen::MatrixXd *H_;

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
	/**
	   TODO:
	   * Initialize the state ekf_.x_ with the first measurement.
	   * Create the covariance matrix.
	   * Remember: you'll need to convert radar from polar to cartesian coordinates.
	   */
	// first measurement
	cout << "EKF: " << endl;
	ekf_.x_ = VectorXd(4);
	ekf_.x_ << 1, 1, 1, 1;

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	    /**
	       Convert radar from polar to cartesian coordinates and initialize state.
	    */
	    float ro;
	    float theta;
	    ro = measurement_pack.raw_measurements_(0);
	    theta = measurement_pack.raw_measurements_(1);
	    ekf_.x_(0) =  ro * std::cos(theta);
	    ekf_.x_(1) =  ro * std::sin(theta);
	    ekf_.Init(ekf_.x_, P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
	    /**
	       Initialize state.
	    */
	    ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
	    ekf_.Init(ekf_.x_, P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
	}
	previous_timestamp_ = measurement_pack.timestamp_;
	// done initializing, no need to predict or update
	is_initialized_ = true;
	return;
    }
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1, 0, dt, 0,
	       0, 1, 0, dt,
	       0, 0, 1, 0,
	       0, 0, 0, 1;
    ekf_.Q_ = MatrixXd(4,4);
    ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
               0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
	       dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
	       0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
    float px = ekf_.x_(0);
    float py = ekf_.x_(1);
    float vx = ekf_.x_(2);
    float vy = ekf_.x_(3);
    float c1 = px*px + py*py;
    float c2 = sqrt(c1);
    float c3 = c1 * c2;
    if(fabs(c1) < 0.0001){
	cout << "CalculateJacobian () - Error - Division by Zero" << endl;
	return;
    } 
    Hj_ << (px/c2), (py/c2), 0, 0,
	    -(py/c1), (px/c1), 0, 0,
	   py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;    

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
       TODO:
       * Update the state transition matrix F according to the new elapsed time.
       - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
       */
    // ekf_set_Time(MeasurementPackage::timestamp);
    cout << "predict begin " << endl;
    ekf_.Predict();
    cout << "predict end " << endl;

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
       TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
       */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	// Radar updates
	ekf_.R_ = R_radar_;
	ekf_.H_ = Hj_;
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
	// Laser updates
	ekf_.R_ = R_laser_;
	ekf_.H_ = H_laser_;
	ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
