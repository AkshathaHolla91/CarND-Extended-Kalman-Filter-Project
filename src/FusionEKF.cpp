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
  
 
  //Initializing the measurement function fro laser
  H_laser_<< 1, 0, 0, 0,
        0, 1, 0, 0;

  

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
   
    
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    //Initializing state transition matrix with an identity matrix, to be updated later
    ekf_.F_=MatrixXd::Identity(4, 4);
    //Creating  the covariance matrix and initializing it
    ekf_.P_=MatrixXd(4,4);
    ekf_.P_<< 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
    

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      
      */
      float rho=measurement_pack.raw_measurements_(0);
      float phi=measurement_pack.raw_measurements_(1);
      float rho_dot=measurement_pack.raw_measurements_(2);
      
      ekf_.x_(0)=rho*cos(phi);
      ekf_.x_(1)=rho*sin(phi);
      ekf_.x_(2)=rho_dot*cos(phi);
      ekf_.x_(3)=rho_dot*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_(0)=measurement_pack.raw_measurements_(0);
      ekf_.x_(1)=measurement_pack.raw_measurements_(1);
      ekf_.x_(2)=0;
      ekf_.x_(3)=0;
      
    }
    previous_timestamp_= measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

 
  //Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  float noise_ax=9;
  float noise_ay=9;
  float delta_t=(measurement_pack.timestamp_-previous_timestamp_)/1000000.0;
  previous_timestamp_=measurement_pack.timestamp_;
  // Update the state transition matrix F according to the new elapsed time.
  ekf_.F_(0,2)=delta_t;
  ekf_.F_(1,3)=delta_t;
  
  float delta_t_square=delta_t* delta_t;
  float delta_t_cube=delta_t_square*delta_t;
  float delta_t_quad=delta_t_cube* delta_t;
  
  //Update the process noise covariance matrix.
  ekf_.Q_=MatrixXd(4,4);
  ekf_.Q_<< delta_t_quad*noise_ax/4, 0 ,delta_t_cube* noise_ax/2, 0,
            0, delta_t_quad*noise_ay/4, 0, delta_t_cube* noise_ay/2,
            delta_t_cube*noise_ax/2, 0, delta_t_square*noise_ax, 0,
            0, delta_t_cube*noise_ay/2, 0, delta_t_square*noise_ay;
  
  

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
  
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tools;
    ekf_.R_=R_radar_;
    ekf_.H_=tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
  ekf_.R_=R_laser_;
  ekf_.H_=H_laser_;
  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
