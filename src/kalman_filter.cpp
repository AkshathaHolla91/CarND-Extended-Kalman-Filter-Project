#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;



KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  
    * predict the state 
  */
  x_=F_ * x_;
  MatrixXd Ft=F_.transpose();
  P_=F_*P_*Ft+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  
    * update the state  by using Kalman Filter equations
  */
  VectorXd y=z - (H_*x_);
  MatrixXd Ht=H_.transpose();
  MatrixXd S_=H_*P_*Ht+R_;
  MatrixXd Si=S_.inverse();
  MatrixXd K=P_*Ht*Si;
  //new estimate
  x_=x_+(K * y);
  long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_=(I-K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
 
    * update the state by using Extended Kalman Filter equations
  */
  float px=x_(0);
  float py=x_(1);
  float vx=x_(2);
  float vy=x_(3);
  //Converting the input from cartesian to polar co-ordinates
  float rho=sqrt(px*px+py*py);
  float phi=atan2(py,px);
	//Check for division by zero and update with small value if zero
  if(rho < 0.00001){
    rho=0.0001;
  }
  float rho_dot=(px*vx+py*vy)/rho;
  VectorXd H_est(3);
  H_est << rho, phi, rho_dot;
  
  VectorXd y=z - H_est;
	//Normalize the angle phi to lie between the range -PI to PI
	float angle=y(1);
	if(angle<-M_PI){
		angle+=2*M_PI;
	}
	else if(angle> M_PI){
		angle-=2*M_PI;
	}
	y(1)=angle;
  MatrixXd Ht=H_.transpose();
  MatrixXd S_=H_*P_*Ht+R_;
  MatrixXd Si=S_.inverse();
  MatrixXd K=P_*Ht*Si;
  
  x_=x_+(K * y);
  long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_=(I-K * H_) * P_;
  
}
