#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define EPS 0.0001

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

  //Finish initializing the FusionEKF.
  noise_ax = 9.;
  noise_ay = 9.;
  


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


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
    cout << "Initialization Beginning" << endl;
    VectorXd x(4);
    //Check first measurement and convert from polar to cartesian if its radar. Store it
    //in x
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      x << rho*cos(phi), rho*sin(phi), 0.f, 0.f;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.f, 0.f;
    }
	previous_timestamp_ = measurement_pack.timestamp_;
	//Initialize initial state covariance matrix. 
	//vx and vy set to 0
    //Use 1000 for covariance to indicate that they are unreliable.
	MatrixXd P(4, 4);
    P << 1, 0, 0, 0,
		 0, 1, 0, 0,
	     0, 0, 1000, 0,
		 0, 0, 0, 1000;
	//Laser measurement matrix initialization
	H_laser_ << 1, 0, 0, 0,
    					 0, 1, 0, 0;
	//Transition
    MatrixXd F(4, 4);
    F << 1, 0, 0, 0,
    	 0, 1, 0, 0,
    	 0, 0, 1, 0,
    	 0, 0, 0, 1;
    	 
    
    
    MatrixXd Q(4,4);
    ekf_.Init( x, /*x_in*/ 
        P, /*P_in*/ 
        F, /*F_in*/
        H_laser_, /*H_in*/ 
        Hj_, /*Hj_in*/ 
        R_laser_, /*R_in*/ 
        R_radar_, /*R_ekf_in*/ 
        Q ); /*Q_in*/
			   
	
    is_initialized_ = true;
        cout << "Initialization Completed" << endl;
    return;
  }

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
   
   //TIME computation, and check to make sure time has elapsed
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  if (dt > 0.0 ){
    // Update the matrix
    ekf_.F_ << 1, 0, dt, 0,
			   0, 1, 0, dt,
			   0, 0, 1, 0,
			   0, 0, 0, 1;
    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;
    float dt4d4 = dt4 / 4.;
    float dt3d2 = dt3 / 2.;
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4d4 * noise_ax, 0, dt3d2 * noise_ax, 0,
	           0, dt4d4 * noise_ay, 0, dt3d2 * noise_ay,
	           dt3d2 * noise_ax, 0, dt2 * noise_ax, 0,
 	           0, dt3d2 * noise_ay, 0, dt2 * noise_ay;
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  //Use Jacobian and call UpdateEKF if updating Radar data, otherwise just call Update
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);
    //ekf_.R_ = R_radar_;
    ekf_.UpdateEKF( measurement_pack.raw_measurements_ );
  } else {
    //ekf_.H_ = H_laser_;
    //ekf_.R_ = R_laser_;
    ekf_.Update( measurement_pack.raw_measurements_ );
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << "Measurement Processed" << endl;
}
