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
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size()==0||estimations.size()!=ground_truth.size()){
        cout<<"CalculateRMSE()-invalid input"<<endl;
        return rmse;
  }
  //accumulate squared residuals
  VectorXd sub(4);
  sub << 0,0,0,0;
  VectorXd sum(4);
  sum << 0,0,0,0;
  for(int i=0; i < estimations.size(); ++i){
	sub=estimations[i]-ground_truth[i];
	sum=sum.array()+sub.array()*sub.array();
  }

  //calculate the mean
  VectorXd mean(4);
  mean << 0,0,0,0;
  mean=sum/estimations.size();
  //calculate the squared root
  rmse=mean.array().sqrt();
  //return the result
  return rmse;  

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float den = sqrt(px*px+py*py);
	

	//check division by zero
	if(fabs(px*px+py*py) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/den), (py/den), 0, 0,
		  -(py/pow(den,2)), (px/pow(den,2)), 0, 0,
		  py*(vx*py - vy*px)/pow(den,3), px*(px*vy - py*vx)/pow(den,3), px/den, py/den;


  return Hj;
}
