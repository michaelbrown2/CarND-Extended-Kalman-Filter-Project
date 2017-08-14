#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#define PI 3.14159265

using namespace std;


using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_in, 
                        MatrixXd &R_ekf_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  R_ = R_in;
  R_ekf_ = R_ekf_in;
  Q_ = Q_in;
  I_ = MatrixXd::Identity( x_.size(), x_.size() );
}


//PREDICT the state 
void KalmanFilter::Predict() {
  //We would need to add +u to the next line if we had motion
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}
//UPDATE the state
void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  Process(y);
}
//UPDATE EKF
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //Avoid dividing by zero:
  if( x_[0] == 0. | x_[1] == 0. )
    return;
  double rho = sqrt( x_(0) * x_(0) + x_(1) * x_(1) );
  double theta = atan2( x_(1), x_(0) );
  double rho_dot = ( x_(0) * x_(2) + x_(1) * x_(3) ) / rho;
  VectorXd h = VectorXd(3);
  h << rho, theta, rho_dot;
  VectorXd y = z - h;
  if( y[1] > PI )
    y[1] -= 2.f*PI;
  if( y[1] < -PI )
    y[1] += 2.f*PI;
  MatrixXd Ht = Hj_.transpose();
  MatrixXd S = Hj_ * P_ * Ht + R_ekf_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  x_ = x_ + ( K * y );
  P_ = ( I_ - K * Hj_ ) * P_;
  
}
//attempt to simplify, both use the same equations (essentially) but I could not
//get the Jacobian to function in the H_ Matrix.  Left it in to further my progress later.
void KalmanFilter::Process(const VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  x_ = x_ + ( K * y );
  P_ = ( I_ - K * H_ ) * P_;
}
