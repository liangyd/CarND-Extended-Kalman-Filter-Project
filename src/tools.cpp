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

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //check division by zero
  if(px==0&&py==0){
    cout<<"CalculateJacobian()-Error-Division by zero"<<endl;
    Hj<<0,0,0,0,
        0,0,0,0,
        0,0,0,0;
  }
  else{
    double den=sqrt(px*px+py*py);
    Hj(0,0)=px/den;
    Hj(0,1)=py/den;
    Hj(1,0)=-py/pow(den,2);
    Hj(1,1)=px/pow(den,2);
    Hj(2,0)=py*(vx*py-vy*px)/pow(den,3);
    Hj(2,1)=px*(vy*px-vx*py)/pow(den,3);
    Hj(2,2)=px/den;
    Hj(2,3)=py/den;
  }
  
  return Hj;
}
