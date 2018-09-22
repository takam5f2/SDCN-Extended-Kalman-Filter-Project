#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
 public:
  
  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix for Normal Kalman filter
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix
   * @param Q_in Process covariance matrix
   */
  void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
	    Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_in, Eigen::MatrixXd &Q_in);

  /**
   * Set Uncertainty matrix of F_.
   * @param dt differential time
   */
  void setUncertainty(const float dt);

  /**
   * Set Process covariance of Q_.
   * @param dt differential time
   * @param noise_ax: noise for ax 
   * @param noise_ay: noise for ay
   */
  void setProcessCov(const float dt, const float noise_ax, const float noise_ay);
  
  /**
   * Set measurement matrix of H_ for EKF.
   */
  void setMeasurement_Radar(void);

  /**
   * Set measurement matrix of H_ for normal KF.
   * @param H_in Measurement matrix for laser
   */
  void setMeasurement_Laser(Eigen::MatrixXd &H_in);

  /**
   * Set measurement matrix of H_ for normal KF.
   * @param R_in Measurement noise matrix
   */
  void setMeasurementNoise(Eigen::MatrixXd &R_in);
  
  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

};

#endif /* KALMAN_FILTER_H_ */
