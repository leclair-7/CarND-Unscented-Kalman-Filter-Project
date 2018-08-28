#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 2.83;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // 
  std_yawdd_ = .2;
  
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
  lambda_ =  3 - n_aug_;

  num_runs_ = 0;

  // starting with Zsig_ uninitialized... 

  /*
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    num_runs_ += 1;

   // cout<< "num_runs_: " << num_runs_<<endl;
    double dt;
    if ( !is_initialized_){

        //compute the time elapsed between the current and previous measurements
        dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
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
           x_(3) =  atan2(vy, vx) ;
           x_(4) = 0.0;
          
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
           x_(2) = 2;
           x_(3) = atan2(x_(1), x_(0)) ;
           x_(4) = 0.0;     
        }

        is_initialized_ = true;
        return;
   }//end initialize if block
   dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
   previous_timestamp_ = meas_package.timestamp_;
   
   Prediction( dt );


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { 
        
        // Clear the prediction
        z_radar_pred_ = VectorXd(3);
        S_radar_pred_ = MatrixXd(3, 3);
        
        // Predict it
        PredictRadarMeasurement( &z_radar_pred_, &S_radar_pred_);

        UpdateRadar( meas_package );

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {     
          UpdateLidar( meas_package ); 
    }
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
  int n_z = 3;

  //define spreading parameter
  double lambda_ = 3 - n_aug_;

  VectorXd weights = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  /*
  Xsig_pred_ 
  x_
  P_
  Zsig_
  z_radar_pred_
  S_radar_pred_ 
  */

  double rho = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double rho_dot = meas_package.raw_measurements_[2];

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z << rho, phi, rho_dot;

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

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
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

  // std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  // std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}


void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda_ = 3 - n_aug_;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  //radar measurement noise standard deviation radius in m
  double std_radr_ = 0.3;

  //radar measurement noise standard deviation angle in rad
  double std_radphi_ = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  double std_radrd_ = 0.1;

  // Matrix with predicted sigma points
  //Xsig_pred_

  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

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
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig_.col(i);
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

    S = S + weights(i) * z_diff * z_diff.transpose();
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
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

    Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    AugmentedSigmaPoints(&Xsig_pred_);
    
    // cout<< "Bubbles:"<<endl;
    // cout<< Xsig_pred_ << endl;
    SigmaPointPrediction(&Xsig_pred_, delta_t);
    
    //cout << "Class var before predict mean x_" << x_ << endl;
    PredictMeanAndCovariance( &x_,  &P_ );
    //cout << "Class var post-predict-mean x_" << x_ << endl;


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
}




void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

  VectorXd weights = VectorXd(2*n_aug_+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);


  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x + weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2. * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2. * M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //std::cout << "Weights state" << std::endl;
  //std::cout << weights << std::endl;

  /*
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;
  */

  //write result
  *x_out = x;
  *P_out = P;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, double delta_t) {

  MatrixXd Xsig_aug = *Xsig_out;

  //cout<< "Bubbles22 (beginning of sigma predict):"<<endl;
  //cout<< Xsig_aug << endl;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  float vk_over_psi_k;
  float psi_k;
  float psi_dot_k;
  float nu_a_k;
  float nu_psi_dotdot_k;
  float vk;
  //2 * n_aug + 1
  for (int i=0; i < 2 * n_aug_ + 1; i++)
  {
     // initialize sigma point for prediction
     MatrixXd predicted_sigma_point = MatrixXd(n_x_,1);
     
     MatrixXd curr_x_k          = MatrixXd(n_x_,1);
     MatrixXd movement_timestep = MatrixXd(n_x_,1);
     MatrixXd process_noise     = MatrixXd(n_x_,1);

     curr_x_k << Xsig_aug(0, i), Xsig_aug(1, i), Xsig_aug(2,i), Xsig_aug(3,i), Xsig_aug(4,i) ;

    vk              = curr_x_k(2);
    psi_k           = curr_x_k(3);
    psi_dot_k       = curr_x_k(4);
    vk_over_psi_k   = vk / psi_dot_k;
    nu_a_k          = Xsig_aug(5, i);
    nu_psi_dotdot_k = Xsig_aug(6, i);

     // if the yaw rate is 0
     if ( Xsig_aug( 4, i ) < .0001){

        movement_timestep << vk * cos( psi_k ) * delta_t ,
                             vk * sin( psi_k ) * delta_t ,
                                  0,
                                  psi_dot_k * delta_t,
                                  0 ;
        
        process_noise << .5 * pow(delta_t, 2) * cos(psi_k) * nu_a_k, 
                         .5 * pow(delta_t, 2) * sin(psi_k) * nu_a_k, 
                          delta_t * nu_a_k,
                          .5 * pow(delta_t, 2) * nu_psi_dotdot_k,
                          delta_t * nu_psi_dotdot_k
                          ;

     } else{

        movement_timestep << vk_over_psi_k*(sin(psi_k+psi_dot_k*delta_t) - sin(psi_k)) ,
                        vk_over_psi_k*(-1 * cos(psi_k+psi_dot_k*delta_t) + cos(psi_k)) ,
                                  0,
                                  psi_dot_k * delta_t,
                                  0 ;
        
        process_noise << .5 * pow(delta_t, 2) * cos(psi_k) * nu_a_k, 
                         .5 * pow(delta_t, 2) * sin(psi_k) * nu_a_k, 
                          delta_t * nu_a_k,
                          .5 * pow(delta_t, 2) * nu_psi_dotdot_k,
                          delta_t * nu_psi_dotdot_k
                          ; 
     }
    
     predicted_sigma_point = curr_x_k + movement_timestep + process_noise;

     Xsig_pred.col(i) = predicted_sigma_point;
  }

  *Xsig_out = Xsig_pred;
  //cout<< "Bubbles (end of sigma predict): "<< delta_t << endl;
  //cout<< Xsig_pred << endl;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {


  MatrixXd Q = MatrixXd(2,2);
  Q << pow(std_a_, 2), 0, 0,pow(std_yawdd_, 2); 


  //create augmented mean vector, put current state values in it, x_ 
  VectorXd x_aug = VectorXd(7);
  for(int i=0; i < n_aug_;i++){
      if (i < n_x_){  
          x_aug(i) = x_(i);
          
      }else{
          x_aug(i)=0;
      }
  }

  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(n_x_,n_x_) = Q(0,0);
  P_aug(n_x_+1,n_x_+1) = Q(1,1);
  //the sqrt of the P_aug matrix  
  MatrixXd A = P_aug.llt().matrixL();
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/
 
  //create augmented mean state
  //create augmented covariance matrix
  //create square root matrix
  //create augmented sigma points
  
  
for(int row=0; row<n_aug_; row++){
    Xsig_aug(row,0) = x_aug(row);
    for(int col=1; col<= n_aug_ ; col++)  {
        Xsig_aug(row,col) = x_aug(row) + sqrt(lambda_ + n_aug_) * A(row,(col-1));
        Xsig_aug(row,col+n_aug_) = x_aug(row) - sqrt(lambda_ + n_aug_) * A(row,(col-1)); 
    }
}
    
  *Xsig_out = Xsig_aug;

}

