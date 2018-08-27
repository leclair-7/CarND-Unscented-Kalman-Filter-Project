#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"
#include "tools.h"

using namespace std;


int main()
{
  

  // Create a Kalman Filter instance
  UKF ukf;

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  
    	//  long long timestamp;
//R 1.014892e+00  5.543292e-01  4.892807e+00  1477010443050000  
        //8.599968e-01  6.000449e-01  5.199747e+00  1.796856e-03  3.455661e-04  1.382155e-02

    	  
    	  string sensor_type = "R";
    	  MeasurementPackage meas_package;

    	  if (sensor_type.compare("L") == 0) {
      	  		meas_package.sensor_type_ = MeasurementPackage::LASER;
          		meas_package.raw_measurements_ = VectorXd(2);
          		//float px;
      	  		//float py;
          		
          		//meas_package.raw_measurements_ << px, py;
          		
          		//meas_package.timestamp_ = timestamp;
          } else if (sensor_type.compare("R") == 0) {

      	  		meas_package.sensor_type_ = MeasurementPackage::RADAR;
          		meas_package.raw_measurements_ = VectorXd(3);
          		float ro = 1.014892e+00;
      	  		float theta = 5.543292e-01;
      	  		float ro_dot = 4.892807e+00;
          		

          		meas_package.raw_measurements_ << ro,theta, ro_dot;
          		

          		meas_package.timestamp_ = 1477010443050000;
          }
        float x_gt = 8.599968e-01;
    	  float y_gt = 6.000449e-01;
    	  float vx_gt = 5.199747e+00;
    	  float vy_gt =  1.796856e-03;
    	  
    	  VectorXd gt_values(4);
    	  gt_values(0) = x_gt;
    	  gt_values(1) = y_gt; 
    	  gt_values(2) = vx_gt;
    	  gt_values(3) = vy_gt;
    	  ground_truth.push_back(gt_values);
          
          //Call ProcessMeasurment(meas_package) for Kalman filter
    	  ukf.ProcessMeasurement(meas_package);    	  


    	  //Push the current estimated x,y positon from the Kalman filter's state vector

        VectorXd estimate(4);

        double p_x = ukf.x_(0);
        double p_y = ukf.x_(1);
        double v  = ukf.x_(2);
        double yaw = ukf.x_(3);

        cout << "Some output:" << endl;
        cout<< p_x << endl;
        cout<< p_y << endl;
        cout<< v << endl;
        cout<< yaw << endl;

/*
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        estimate(0) = p_x;
        estimate(1) = p_y;
        estimate(2) = v1;
        estimate(3) = v2;
        
        estimations.push_back(estimate);

        VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);
*/


//R 1.047505e+00  3.892401e-01  4.511325e+00  1477010443150000 
// 1.379955e+00  6.006288e-01  5.198979e+00  1.077814e-02  2.073124e-03  2.763437e-02

        meas_package.sensor_type_ = MeasurementPackage::RADAR;
        meas_package.raw_measurements_ = VectorXd(3);
        float ro2 = 1.047505e+00 ;
        float theta2 = 3.892401e-01;
        float ro_dot2 = 4.511325e+00;
        

        meas_package.raw_measurements_ << ro2,theta2, ro_dot2;
        

        meas_package.timestamp_ = 1477010443150000;
          
        float x_gt2 = 1.379955e+00;
        float y_gt2 = 6.006288e-01;
        float vx_gt2 = 5.198979e+00;
        float vy_gt2 =  1.077814e-02;
        
        //VectorXd gt_values(4);
        gt_values(0) = x_gt2;
        gt_values(1) = y_gt2; 
        gt_values(2) = vx_gt2;
        gt_values(3) = vy_gt2;
        ground_truth.push_back(gt_values);
          
          //Call ProcessMeasurment(meas_package) for Kalman filter
        ukf.ProcessMeasurement(meas_package);


        double p_x2 = ukf.x_(0);
        double p_y2 = ukf.x_(1);
        double v2  = ukf.x_(2);
        double yaw2 = ukf.x_(3);

        cout << "Some output:" << endl;
        cout<< p_x2 << endl;
        cout<< p_y2 << endl;
        cout<< v2 << endl;
        cout<< yaw2 << endl;


}
  





















































































