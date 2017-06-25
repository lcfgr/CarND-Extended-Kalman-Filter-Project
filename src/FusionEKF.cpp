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

  ekf_.x_ = VectorXd(4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.Q_  = MatrixXd(4,4);
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);



  ekf_.x_ << 1, 1, 1, 1;

  //the initial state covariance matrix P_
  ekf_.P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  //the initial state transition matrix F_
  ekf_.F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  //the initial laser measurement matrix H_laser_
  H_laser_ << 1, 0, 0, 0,
          0, 1, 0, 0;

  //the initial laser measurement covariance matrix R_laser_
  R_laser_ << 0.0225, 0,
          0, 0.0225;

  //the initial radar measurement covariance matrix R_radar_
  R_radar_ << 0.09, 0, 0,
          0, 0.0009, 0,
          0, 0, 0.09;

  //the initial process covariance matrix Q_
  ekf_.Q_ << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;

  //set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;


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
    cout << "EKF: " << endl;
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      float rho_ = measurement_pack.raw_measurements_[0];
      float phi_ = measurement_pack.raw_measurements_[1];
      float rho_dot_ = measurement_pack.raw_measurements_[2];

      ekf_.x_ <<  rho_ * cos(phi_), rho_ * sin(phi_), rho_dot_ * cos(phi_), rho_dot_ * sin(phi_);


    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      // Don't take null measurements into account
      if (measurement_pack.raw_measurements_[0]  ==0||measurement_pack.raw_measurements_[1]==0){
        cout << "Null laser measurements"<<endl;
        return;
      }

      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // advance timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt * dt_2;
  float dt_4 = dt * dt_3;

  // update F state transition matrix with the time between states
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  //update Q process covariance matrix
  ekf_.Q_ <<  (dt_4/4.0)*noise_ax, 0, (dt_3/2.0)*noise_ax, 0,
          0, (dt_4/4.0)*noise_ay, 0, (dt_3/2.0)*noise_ay,
          (dt_3/2.0)*noise_ax, 0, dt_2*noise_ax, 0,
          0, (dt_3/2.0)*noise_ay, 0, dt_2*noise_ay;

  // call Kalman filter prediction function
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

    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates

    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
