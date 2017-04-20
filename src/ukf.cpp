#include "ukf.h"
#include "tools.h"
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

  n_x_ = 5;

  n_aug_ = 7;

  n_sigma_ = 2 * n_aug_ + 1;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(n_sigma_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2 * n_aug_) = MatrixXd::Constant(2 * n_aug_, 1, 1 / (2 * (lambda_ + n_aug_)));

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  Xsig_pred_.fill(0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; // 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; // 30;

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
}

UKF::~UKF() {}


void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
  if (!is_initialized_) {
    Initialize(measurement_package);
    return;
  }

  double delta_t = (measurement_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = measurement_package.timestamp_;

  this->Predict(delta_t);
  if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
    this->UpdateRadar(measurement_package);
  } else {
    this->UpdateLidar(measurement_package);
  }

}


void UKF::Predict(double delta_t) {
  MatrixXd augmented_sigma_points = this->GenerateSigmaPoints();
  this->PredictSigmaPoints(augmented_sigma_points, delta_t);
  this->PredictNewState();
}

void UKF::PredictNewState() {
  // Predict the state first.
  VectorXd new_mean(n_x_);
  new_mean.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    new_mean += weights_(i) * Xsig_pred_.col(i);
  }
  x_ << new_mean;

  // Then we can predict covariance basing on predicted state.
  MatrixXd X_sub_mean = VectorXd(Xsig_pred_.rows());
  MatrixXd new_covariance(n_x_, n_x_);
  new_covariance.fill(0);
  for (int i=0; i < Xsig_pred_.cols(); i++) {
    X_sub_mean = Xsig_pred_.col(i) - x_;

    while (X_sub_mean(3)> M_PI) X_sub_mean(3)-=2.*M_PI;
    while (X_sub_mean(3)< -M_PI) X_sub_mean(3)+=2.*M_PI;

    new_covariance += weights_(i) * (X_sub_mean * X_sub_mean.transpose());
  }
  P_ << new_covariance;
}

void UKF::PredictSigmaPoints(const MatrixXd &augmented_sigma_points, double delta_t) {
  double delta_t_sq = delta_t * delta_t;
  VectorXd sigma_point(n_aug_);
  VectorXd noise(n_x_);
  VectorXd sigma_point_prediction(n_x_);

  for (int i = 0; i < augmented_sigma_points.cols(); i++) {
    sigma_point << augmented_sigma_points.col(i);
    if (fabs(sigma_point(4)) < 0.001) {
      sigma_point_prediction <<
                             sigma_point(2) * cos(sigma_point(3)) * delta_t,
                             sigma_point(2) * sin(sigma_point(3)) * delta_t,
                             0,
                             0, // yawd eq 0, so entire term is also 0
                             0;
    } else {
      sigma_point_prediction <<
                             (sigma_point(2) / sigma_point(4)) *
                             (sin(sigma_point(3) + sigma_point(4) * delta_t) - sin(sigma_point(3))),
                             (sigma_point(2) / sigma_point(4)) *
                             (-cos(sigma_point(3) + sigma_point(4) * delta_t) + cos(sigma_point(3))),
                             0,
                             sigma_point(4) * delta_t,
                             0;
    }
    noise <<
          0.5 * delta_t_sq * cos(sigma_point(3)) * sigma_point(5) ,
          0.5 * delta_t_sq * sin(sigma_point(3)) * sigma_point(5),
          delta_t * sigma_point(5),
          0.5 * delta_t_sq * sigma_point(6),
          delta_t * sigma_point(6);
    Xsig_pred_.col(i) << sigma_point.head(n_x_) + sigma_point_prediction + noise;
  }
}


MatrixXd UKF::GenerateSigmaPoints() {
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Augmented mean state
  x_aug.head(x_.size()) = x_;
  x_aug.tail(2) << 0, 0;

  // Augmented cov matrix
  P_aug.fill(0);
  P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
  P_aug.bottomRightCorner(2, 2) << std_a_ * std_a_, 0                    ,
                                   0              , std_yawdd_ * std_yawdd_;

  double sqrt_lambda_nx = sqrt(lambda_ + n_aug_);
  MatrixXd sqrtP = P_aug.llt().matrixL();

  Xsig_aug.fill(0.0);
  Xsig_aug.col(0) << x_aug;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) << x_aug + sqrt_lambda_nx * sqrtP.col(i);
    Xsig_aug.col(n_aug_ + i + 1) << x_aug - sqrt_lambda_nx * sqrtP.col(i);
  }

  return Xsig_aug;
}


void UKF::Initialize(const MeasurementPackage &measurement_package) {
  VectorXd initial_x(5);
  if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
      double px, py;
      tie(px, py) = Tools::PolarToCartesian(
              measurement_package.raw_measurements_[0], measurement_package.raw_measurements_[1]
      );
      initial_x << px, py, 0, 0, 0;
    }
    else if (measurement_package.sensor_type_ == MeasurementPackage::LASER) {
      initial_x << measurement_package.raw_measurements_[0], measurement_package.raw_measurements_[1], 0, 0, 0;
    }
  x_ = initial_x;
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 100, 0, 0,
        0, 0, 0, 100, 0,
        0, 0, 0, 0, 100;

  this->time_us_ = measurement_package.timestamp_;
  this->is_initialized_ = true;
}


void UKF::UpdateLidar(MeasurementPackage measurement_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

void UKF::UpdateRadar(MeasurementPackage measurement_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
