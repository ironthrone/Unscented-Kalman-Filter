#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  // Parameters above this line are scaffolding, do not modify

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  n_sigma_ = 2 * n_aug_ + 1;
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  Xsig_pred_.fill(0);

  Q_ = MatrixXd(2, 2);
  Q_ << std_a_ * std_a_, 0,
      0, std_yawdd_ * std_yawdd_;

  R_lidar = MatrixXd(2, 2);
  R_lidar << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  R_radar = MatrixXd(3, 3);
  R_radar << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

  weights_ = VectorXd(n_sigma_);
  weights_(0) = lambda_ * 1.0 / (lambda_ + n_aug_);
  for (int i = 1; i < n_sigma_; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  lidar_nis_file_.open("../lidar_nis_file.txt", ios::out | ios::trunc);
  radar_nis_file.open("../radar_nis_file.txt", ios::out | ios::trunc);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "init by lidar" << endl;
      x_ << meas_package.raw_measurements_(0),
          meas_package.raw_measurements_(1),
          0, 0, 0;
      P_ << 0.1, 0, 0, 0, 0,
          0, 0.1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    } else {
      cout << "init by radar" << endl;
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);
      x_ << rho * cos(phi),
          rho * sin(phi),
          0,
          0,
          0;
      P_ << 0.1, 0, 0, 0, 0,
          0, 0.1, 0, 0, 0,
          0, 0, 10, 0, 0,
          0, 0, 0, 0.1, 0,
          0, 0, 0, 0, 1;

    }
    cout << "init radar x:\n" << x_ << endl;
    is_initialized_ = true;
    last_timestamp = meas_package.timestamp_;
    return;
  }

  double delta = (meas_package.timestamp_ - last_timestamp) / 1000000.0;
  cout << "delta : " << delta << endl;
  last_timestamp = meas_package.timestamp_;
  Prediction(delta);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0);

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(5) = x_;
  cout << "x aug: " << x_aug << endl;
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;

  Xsig_aug.col(0) = x_aug;
  MatrixXd A_aug = P_aug.llt().matrixL();
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(3) * A_aug.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - sqrt(3) * A_aug.col(i);
  }

  for (int i = 0; i < n_sigma_; i++) {
    VectorXd xsig_aug = Xsig_aug.col(i);
    double v_k = xsig_aug(2);
    double yaw_k = xsig_aug(3);
    double yaw_dk = xsig_aug(4);
    double nu_ak = xsig_aug(5);
    double nu_yawddk = xsig_aug(6);

    VectorXd x_delta(5);
    VectorXd x_noise(5);
    // caution use fabs() to check if near 0,or last half will be wired
    if (fabs(yaw_dk) < 0.00001) {
      x_delta << v_k * delta_t * cos(yaw_k), v_k * delta_t * sin(yaw_k), 0, 0, 0;
    } else {
      x_delta << v_k / yaw_dk * (sin(yaw_k + yaw_dk * delta_t) - sin(yaw_k)),
          v_k / yaw_dk * (-cos(yaw_k + yaw_dk * delta_t) + cos(yaw_k)),
          0,
          yaw_dk * delta_t,
          0;
    }
    x_noise << 0.5 * nu_ak * delta_t * delta_t * cos(yaw_k),
        0.5 * nu_ak * delta_t * delta_t * sin(yaw_k),
        nu_ak * delta_t,
        0.5 * nu_yawddk * delta_t * delta_t,
        nu_yawddk * delta_t;
    Xsig_pred_.col(i) = xsig_aug.head(5) + x_delta + x_noise;
  }


  x_ = (Xsig_pred_ * weights_).rowwise().sum();
  P_.fill(0);
  for (int i = 0; i < n_sigma_; i++) {
    VectorXd delta = Xsig_pred_.col(i) - x_;
    //TODO angle normal , why operation to angle above do not need angle normalization
    delta(3) = atan2(sin(delta(3)), cos(delta(3)));
    P_ += weights_(i) * delta * delta.transpose();
  }

//  cout << "predict x : " << x_ << endl;
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

  MatrixXd H(2, 5);
  H << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;
  VectorXd y = meas_package.raw_measurements_ - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R_lidar;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(5, 5);
  P_ = (I - K * H) * P_;

  double epsilon = y.transpose() * S.inverse() * y;
  lidar_nis_file_ << epsilon << endl;
//  cout << "update laser x : " << x_ << endl;
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
  MatrixXd Zsig = MatrixXd(3, n_sigma_);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(3);
  z_pred.fill(0);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(3, 3);
  S.fill(0);

  for (int i = 0; i < n_sigma_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double yaw_d = Xsig_pred_(4, i);

    double rho = hypot(px, py);
    double phi = atan2(py, px);
    double rho_d;
    if (rho < 0.00001) {
      rho_d = 0;
    } else {
      rho_d = (px * v * cos(yaw) + py * v * sin(yaw)) / rho;
    }

    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2, i) = rho_d;

    z_pred += weights_(i) * Zsig.col(i);
  }

  for (int i = 0; i < n_sigma_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R_radar;

  cout << "S for radar: " << S << endl;

  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0);
  for (int i = 0; i < n_sigma_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // different with update laser
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // different with update laser
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K(n_x_, 3);
  K = Tc * S.inverse();
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  // different with update laser
  z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)));
  x_ = x_ + K * z_diff;
  MatrixXd i;
  i = MatrixXd::Identity(n_x_, n_x_);
  P_ = P_ - K * S * K.transpose();

  double epsilon = z_diff.transpose() * S.inverse() * z_diff;
  radar_nis_file << epsilon << endl;
//  cout << "update radar x : " << x_ << endl;
}
