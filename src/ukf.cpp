#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPSILON .001
/**
 * Initializes Unscented Kalman filter
 * 
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // 2 * sqrt 2, because 8 m2/s4 is super fast and maybe humanly possible
  std_a_ = .9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .7;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  is_initialized_ = false;

  previous_timestamp_ = 0;

  n_x_    = 5;
  n_aug_  = 7;
  lambda_ =  5 - n_aug_;

  num_runs_ = 0;

  nis_ = 0.0;

  H_laser_ = MatrixXd(2,5);
  H_laser_.fill(0);
  H_laser_(0,0) = 1;
  H_laser_(1,1) = 1;

  weights_ = VectorXd(2*n_aug_+1);

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  R_ = MatrixXd(2,2);
  R_ <<    std_laspx_*std_laspx_, 0, 
          0, std_laspy_*std_laspy_;

}

UKF::~UKF() {}

double UKF::Safe_atan2(double vy, double vx){
  if (vy < EPSILON and vx < EPSILON){
    return 0.0;
  } else {
    return atan2(vy, vx);
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::InitValues(MeasurementPackage meas_package){

        //compute the time elapsed between the current and previous measurements
        previous_timestamp_ = meas_package.timestamp_;

        x_ = VectorXd(n_x_);
        P_ = MatrixXd::Identity(n_x_,n_x_);


        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {     

           
           double rho = meas_package.raw_measurements_[0];
           double phi = meas_package.raw_measurements_[1];
           double rho_dot = meas_package.raw_measurements_[2];
           double x= rho * cos(phi);
           double y= rho * sin(phi);
           double vx = rho_dot * cos(phi);
           double vy = rho_dot * sin(phi);
           
           double velocity = sqrt( pow(vx,2) + pow(vy,2) );
           
           /*
            px
            py
            v
            yaw
            yaw_rate_of_change
           */

           x_(0) = x;
           x_(1) = y;
           x_(2) = velocity;
           x_(3) = Safe_atan2(vy, vx);     
           x_(4) = 0.2;
          
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {     

           /*
            px
            py
            v
            yaw
            yaw_rate_of_change
           */

           x_(0) = meas_package.raw_measurements_(0);
           x_(1) = meas_package.raw_measurements_(1);
           x_(2) = 1.5;           
           x_(3) = Safe_atan2(x_(1), x_(0)) ;         
           x_(4) = 0.0;  

        }
      return;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    nis_array_[num_runs_] = nis_;
    /*
    num_runs_ += 1;
    cout<< num_runs_<<endl;
    if (num_runs_ == 499)
    {
        ofstream loggerman("nis_array_.log");
        
        if (loggerman.is_open())
        {
          for(int count = 0; count < num_runs_; count ++){
              loggerman << nis_array_[count] << "\n" ;
          } 
        }
        loggerman.close();
    }
    */
    //cout<< "num_runs_: " << num_runs_<<endl;
    if ( !is_initialized_){
        InitValues(meas_package);

        is_initialized_ = true;
        return;
   }
   
   double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
   previous_timestamp_ = meas_package.timestamp_;
   
   Prediction( dt );


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_) { 
        
        // Clear the prediction
        z_radar_pred_ = VectorXd(3);
        S_radar_pred_ = MatrixXd(3, 3);
        
        // Predict it
        PredictRadarMeasurement( &z_radar_pred_, &S_radar_pred_);

        UpdateRadar( meas_package, z_radar_pred_, S_radar_pred_ );

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_) {     
      
       UpdateLidar( meas_package );
      }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    AugmentedSigmaPoints(&Xsig_pred_);
    
    SigmaPointPrediction(&Xsig_pred_, delta_t);
    
    PredictMeanAndCovariance( &x_,  &P_ );
}


void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t) {


  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred.fill(0.0);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {    
    double p_x      = Xsig_pred_(0,i);
    double p_y      = Xsig_pred_(1,i);
    double v        = Xsig_pred_(2,i);
    double yaw      = Xsig_pred_(3,i);
    double yawd     = Xsig_pred_(4,i);
    double nu_a     = Xsig_pred_(5,i);
    double nu_yawdd = Xsig_pred_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  *Xsig_out = Xsig_pred;
}
/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package, VectorXd z_radar_pred_, MatrixXd S_radar_pred_) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  // already have sig pred, i.e. Xsig_pred_

  /*
  Xsig_pred_ 
  x_
  P_
  Zsig_
  z_radar_pred_
  S_radar_pred_ 
  */

  double rho     = meas_package.raw_measurements_[0];
  double phi     = meas_package.raw_measurements_[1];
  double rho_dot = meas_package.raw_measurements_[2];

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z(0) = rho;
  z(1) = phi;
  z(2) = rho_dot;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_radar_pred_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_radar_pred_.inverse();

  //residual
  VectorXd z_diff = z - z_radar_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_radar_pred_ * K.transpose();

  nis_ = (z - z_radar_pred_).transpose() * ( S_radar_pred_.inverse()) * (z - z_radar_pred_);
  //std::cout<< "radar nis_: " << nis_ << endl;

   //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
   //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
}


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = Safe_atan2(p_y,p_x);                         //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / std::max(EPSILON, sqrt(p_x*p_x + p_y*p_y));   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig_.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred;

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

  
/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  /*
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;
  */

  //write result
  *z_out = z_pred;
  *S_out = S;
}


void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
  /*
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;
  std::cout <<" Xsig_pred_"<<std::endl;
  std::cout << Xsig_pred_ << std::endl;
  */

  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);


  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }
  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2. * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2. * M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //std::cout << "Weights state" << std::endl;
  //std::cout << weights_ << std::endl;

  //std::cout << "Weights M_PI" << std::endl;
  //std::cout << M_PI << std::endl;

  
  /*
  */

  //write result
  *x_out = x;
  *P_out = P;
}


void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {


  MatrixXd Q = MatrixXd(2,2);
  Q << pow(std_a_, 2), 0, 0,pow(std_yawdd_, 2); 


  //create augmented mean vector, put current state values in it, x_ 
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //populate x_aug
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //populate P_aug
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = Q(0, 0);
  P_aug(n_x_ + 1, n_x_ + 1) = Q(1, 1);
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
    
  *Xsig_out = Xsig_aug;
}

void UKF::PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
  /*
    Since Lidar is linear, we don't need the sigma points, given nonlinear data that would be relevant
  */

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
  /**
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Uses lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  Calculation of NIS done in ProcessMeasurement
  */
  int n_z = 2;

  double px     = meas_package.raw_measurements_[0];
  double py     = meas_package.raw_measurements_[1];
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z(0) = px;
  z(1) = py;
  
  //1
  VectorXd y_ = (z - (H_laser_ * x_));
  //2
  MatrixXd S_ = (H_laser_ * P_ * (H_laser_.transpose())) + R_;
  //3
  MatrixXd K_ = (P_ * (H_laser_.transpose())  * (S_.inverse())) ;
  //4
  x_ = x_ + K_ * y_;
  //5
  MatrixXd I = MatrixXd::Identity(5,5);
  //6
  P_ = (I - K_ * H_laser_ ) * P_;
 
}