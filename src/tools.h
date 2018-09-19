#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
  * A helper method to calculate Process covariance.
  */
  MatrixXd CalculatePCovariance(const float delta_t, const float noise_ax, const float noise_ay);

  /**
  * A helper method to Prediction matrix
  */
  MatrixXd PredictionMatrix(const float delta_t);
  

};

#endif /* TOOLS_H_ */
