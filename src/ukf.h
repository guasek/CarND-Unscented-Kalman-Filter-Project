#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;

    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;

    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;

    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;

    ///* state covariance matrix
    MatrixXd P_;

    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;

    ///* time when the state is true, in us
    long long time_us_;

    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;

    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;

    ///* Laser measurement noise standard deviation position1 in m
    double std_laspx_;

    ///* Laser measurement noise standard deviation position2 in m
    double std_laspy_;

    ///* Radar measurement noise standard deviation radius in m
    double std_radr_;

    ///* Radar measurement noise standard deviation angle in rad
    double std_radphi_;

    ///* Radar measurement noise standard deviation radius change in m/s
    double std_radrd_;

    ///* Weights of sigma points
    VectorXd weights_;

    ///* State dimension
    int n_x_;

    ///* Augmented state dimension
    int n_aug_;

    ///* Number of sigma points
    int n_sigma_;

    int radar_z_dim_;

    int lidar_z_dim;

    ///* Sigma point spreading parameter
    double lambda_;

    ///* the current NIS for radar
    double NIS_radar_;

    ///* the current NIS for laser
    double NIS_laser_;

    /**
     * Constructor
     */
    UKF();

    /**
     * Destructor
     */
    virtual ~UKF();

    /**
     * ProcessMeasurement
     *
     * @param MeasurementPackage meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);

    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     *
     * @param double delta_t Time between k and k+1 in s
     */
    void Predict(double delta_t);

    /**
     * Updates the state and the state covariance matrix using a laser measurement
     *
     * @param MeasurementPackage measurement_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage measurement_package);

    /**
     * Updates the state and the state covariance matrix using a radar measurement
     *
     * @param MeasurementPackage measurement_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage measurement_package);

    /**
     * Performs KF initialization.
     *
     * @param MeasurementPackage& measurement_package Measurement package.
     */
    void Initialize(const MeasurementPackage &measurement_package);

    /**
     * Generates sigma point matrix.
     *
     * @param Xsig_out
     */
    MatrixXd GenerateSigmaPoints();

    /**
     * Does a sigma point prediction basing on augmented state vector.
     *
     * @param MatrixXd& augmented_sigma_points
     */
    void PredictSigmaPoints(const MatrixXd &augmented_sigma_points, double delta_t);

    /**
     * Does mean and covariance prediction basing on sigma point predictions.
     */
    void PredictNewState();

    /**
     * Updates mean vector and covariance matrix.
     */
    void UpdateMeanAndCovariance(
            const MatrixXd &R,
            const MatrixXd &Zsig,
            int measurement_z_dim,
            MeasurementPackage &measurement_package
    );
};
#endif /* UKF_H */
