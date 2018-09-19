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

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
        0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    // process covariance matrix
    H_laser_ << 1, 0, 0, 0,
	        0, 1, 0, 0;
    
    // Uncertainty covariance matrix
    P_ = MatrixXd(4,4);
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
	    x_state_(0) =  ro * std::cos(theta);
	    x_state_(1) =  ro * std::sin(theta);
	    ekf_.Init(x_state_, P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
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
    Tools tools;
    // Hj_ = tools.CalculationJacobian(x_state)

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
