#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
    return rmse;
  }

  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    rmse = rmse + (residual.array() * residual.array()).matrix();
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

std::tuple<double, double> Tools::PolarToCartesian(double rho, double phi) {

  double px = rho * cos(phi);
  double py = rho * sin(phi);

  return std::make_tuple(px, py);
};

double Tools::NormalizeAngle(double angle) {
  while (angle > M_PI) {
    angle -= 2.*M_PI;
  }
  while (angle < -M_PI) {
    angle += 2.*M_PI;
  }
  return angle;
}