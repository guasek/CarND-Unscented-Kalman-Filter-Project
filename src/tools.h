#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include <iostream>

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * Converts polar coordinates to cartesian.
   */
  static std::tuple<double, double> PolarToCartesian(double rho, double phi);

  /**
   * Normalizes given angle so it's between -pi,pi.
   */
  static double NormalizeAngle(double angle);
};

#endif /* TOOLS_H_ */
