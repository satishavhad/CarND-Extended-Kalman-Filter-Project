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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
	          0, 1, 0, 0;

  Hj_ <<	0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0;
  // Defind initial covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1000, 0,
	  0, 0, 0, 1000;
  
  //Initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
	  0, 1, 0, 1,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

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
	// Defind initial covariance matrix
	  ekf_.P_ = MatrixXd(4, 4);
	  ekf_.P_ <<  1, 0, 0, 0,
				  0, 1, 0, 0,
				  0, 0, 1000, 0,
				  0, 0, 0, 1000;
	
	
    ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		// Initialize the Radar sensor input
		float rho = measurement_pack.raw_measurements_[0];
		float phi = measurement_pack.raw_measurements_[1];
		float rhodot = measurement_pack.raw_measurements_[2];
		
		// Covert radar from polar to cartesian coordinates
		float px = rho * cos(phi);
		float py = rho * sin(phi);
		float vx = rhodot * cos(phi);
		float vy = rhodot * sin(phi);
		ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		//Velocity is not captured by Laser sensors
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;
    }

	// first measurement
	//cout << "EKF: " << ekf_.x_ << endl;

	// Capture timestamp
	previous_timestamp_ = measurement_pack.timestamp_;
	
	//wait till px reading is not zero
	if (measurement_pack.raw_measurements_[0] == 0) {
		return;
	}

    // done initializing, no need to predict or update
    is_initialized_ = true;
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

  //compute delta time between measurements
  float deltaT = (measurement_pack.timestamp_ - previous_timestamp_);
  // Time is measured in seconds
  deltaT = deltaT / 1000000.0; 
  
  // Capture timestamp
  previous_timestamp_ = measurement_pack.timestamp_;

  // Update the state transition matrix F according to the new elapsed time
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, deltaT, 0,
			 0, 1, 0, deltaT,
			 0, 0, 1, 0,
			 0, 0, 0, 1;
  
  // Update the process noise covariance matrix	
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  
  float dt_sqar = deltaT * deltaT; //DeltaT Square
  float dt_cube = deltaT * deltaT * deltaT; //DeltaT Cube
  float delta4 = deltaT * deltaT * deltaT * deltaT; //DeltaT raise to 4
  float delta4_4 = delta4 / 4;
  float dt_cube_2 = dt_cube / 2;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<	delta4_4 * noise_ax,	0,					dt_cube_2 * noise_ax,	0,
				0,						delta4_4 * noise_ay,	0,					dt_cube_2 * noise_ay,
				dt_cube_2 * noise_ax,	0,					dt_sqar * noise_ax,		0,
				0,						dt_cube_2 * noise_ay,	0,					dt_sqar * noise_ay;


  ekf_.Predict();

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
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
