#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

// predict
void KalmanFilter::Predict() {
  x_ = F_ * x_ ;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

//upadte the Kalman filter based on the Laser measurement
void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  // New state
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
     update the state by using Extended Kalman Filter equations
  */

  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float phi = atan2(x_(1), x_(0));
  float rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  // normalize the phi value in the y vector.
  while ( y(1) < -M_PI || y(1) > M_PI ) {
    if ( y(1) > M_PI ) {
      y(1) -= M_PI;
    } else {
      y(1) += M_PI;
    }
  }
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  // New state
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
