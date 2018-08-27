#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0.0,0.0,0.0,0.0;

	if (estimations.size() != ground_truth.size() || estimations.size() == 0 )
	{
		cout<< "Invalid estimations of ground_truth size"<<endl;
		return rmse;
	}
			
	VectorXd temp(4);

	for(unsigned int i=0; i < estimations.size(); ++i){
		temp = (estimations[i] - ground_truth[i] );
		rmse = rmse.array() + (temp.array() * temp.array() );
	}
	
	rmse = (1.0/estimations.size()) * rmse;
	rmse = rmse.array().sqrt();

	return rmse;


}