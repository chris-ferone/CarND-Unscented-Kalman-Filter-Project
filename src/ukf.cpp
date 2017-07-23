#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  

  
  //* State dimension
  n_x_ = 5;

  //* Augmented state dimension
  n_aug_ = 7;
  
  		//create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  //* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  
  //* Weights of sigma points
  
  weights_ = VectorXd(2*n_aug_+1);
  
  double weight_0 = lambda_/(lambda_+n_aug_);
  //cout << "34563246" << endl;
  weights_(0) = weight_0;
  //cout << "tyjgyj" << endl; 
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
     
  previous_timestamp_ = 0;
  
  dt = 0;
 
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (previous_timestamp_ == 0)
	{
		dt = 0.05;
	}
	else
	{
  dt =  (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
	}
  
   if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	  // cout << "before radar prediction" << endl;
	   //Initalize x_
	    float rho = meas_package.raw_measurements_(0);
		float phi = meas_package.raw_measurements_(1);
		float rhodot = meas_package.raw_measurements_(2);
		x_(0) = rho / sqrt(1 + pow(tan(phi),2));
		x_(1) = x_(0) * tan(phi);
	   
	   Prediction(dt);
	   //cout << "before radar update" << endl;
	   UpdateRadar(meas_package);
	   //cout << "after radar update" << endl;
	   
   }
   else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	   //cout << "before lidar prediction" << endl;
	   //Initalize x_
	   x_(0) = meas_package.raw_measurements_(0);
	   x_(1) = meas_package.raw_measurements_(1);
	   
	   Prediction(dt);
	   //cout << "before lidar update" << endl;
	   UpdateLidar(meas_package);
	   //cout << "after lidar update" << endl;
	   
   }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  ////////////////////////////////////
  //    Create Sigma Points  /////////
  ////////////////////////////////////
  
    
  
    
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
 //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  //cout << "sigma points created" << endl;
  
  /////////////////////////////////////
  //    Predict Sigma Points  /////////
  /////////////////////////////////////
  
   //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
	//cout << "1" << endl;
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
	//cout << "2" << endl;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
	
	//cout << "3" << endl;
	
	
    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
	//cout << "4" << endl;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
	//cout << "5" << endl;
	//cout << "i:" << i << endl;
	//cout << "px_p" << px_p << endl;
	
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
	//cout << "6" << endl;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
	
  }
	//cout << "sigma points predicted" << endl;
	
	//////////////////////////////////////////////////////
	////     Predicted State x and Covaraiance P     /////
	//////////////////////////////////////////////////////
	
	//create vector for weights
  //VectorXd weights = VectorXd(2*n_aug_+1);
  
  //create vector for predicted state
  //VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  //MatrixXd P = MatrixXd(n_x, n_x);
	
	  // set weights
 

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
	//cout << " state mean predicted" << endl;
	
  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
	
	//cout << " state covariance predicted" << endl;
    //cout << M_PI << endl;
	// state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
	//cout << " before angle normalization" << endl;
    //angle normalization
    while (x_diff(3)> M_PI) {
		//cout << " angle too large" << endl;
		//cout << x_diff(3) << endl;
		x_diff(3)-=2.*M_PI;
	}
    while (x_diff(3)<-M_PI) {
		//cout << " angle too small" << endl;
		//cout << x_diff(3) << endl;
		x_diff(3)+=2.*M_PI;
	}
	//cout << " after angle normalization" << endl;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  //cout << "x and P predicted" << endl;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
    ///////////////////////////////////////////////////////////////////
  //////////////// Predict Lidar Measurement ///////////////////////
  /////////////////////////////////////////////////////////////
  
  //transform sigma points into measurement space
    
	// measurement model
    MatrixXd H_laser_ = MatrixXd(2, 5);
	//measurement matrix
	H_laser_ << 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0;
				
	//measurement covariance matrix - laser
	MatrixXd R_laser_ = MatrixXd(2, 2);
	R_laser_ << 0.0225, 0,
				0, 0.0225;
	
	VectorXd z = VectorXd(2);		
	z(0) = meas_package.raw_measurements_(0);
	z(1) = meas_package.raw_measurements_(1);
				
	//cout << " 1" << endl;
	VectorXd y = z - H_laser_ * x_;
	//cout << " 2" << endl;
	MatrixXd Ht = H_laser_.transpose();
	//cout << " 3" << endl;
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	//cout << " 4" << endl;
	MatrixXd Si = S.inverse();
	//cout << " 5" << endl;
	MatrixXd K =  P_ * Ht * Si;
	//cout << " 6" << endl;
	x_ = x_ + (K * y);
	//cout << " 7" << endl;
	MatrixXd I;
	//cout << " 8" << endl;
	I = MatrixXd::Identity(5, 5);
	//cout << " 9" << endl;
	//cout << sizeof(I) / sizeof(I[0]) << endl;
	//cout << sizeof(K) << endl;
	//cout << sizeof(H_laser_) << endl;
	//cout << sizeof(P_) << endl;
	P_ = (I - K * H_laser_) * P_;
	//cout << " 10" << endl;
  
	  cout << "x: " << x_ << endl;
  cout << "P: " << P_ << endl;
  
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  ///////////////////////////////////////////////////////////////////
  //////////////// Predict Radar Measurement ///////////////////////
  /////////////////////////////////////////////////////////////
  
  //create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  int n_z = 3;
  
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  
  //set vector for weights
  /* VectorXd weights = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  } */
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  
  VectorXd z_diff;
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  ///////////////////////////////////////////////////////////
  ///////////   	UKF update 			/////////////////////
  ///////////////////////////////////////////////////////////
  
   //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
    //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    /* //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI; */

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
	
  VectorXd z = meas_package.raw_measurements_;
	
	
  //residual
  //VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  cout << "x: " << x_ << endl;
  cout << "P: " << P_ << endl;
  
}
