#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double NEAR_ZERO = 0.0001;

double ConvertToPiRange(double angle) {
  while ((angle < -1. * M_PI) || (angle > M_PI)) {
     if (angle < -1. * M_PI) {
       angle += 2. * M_PI;
     } else {
       angle -= 2. * M_PI;
     }
   }
  return angle;
}

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
  std_yawdd_ = 0.8;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x_ = 5;
  //set augmented dimension
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  time_us_ = 0; // Timestamp initialization
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5F / (lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/

  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double range = meas_package.raw_measurements_[0];
      double bearing = meas_package.raw_measurements_[1];
      double px = range * cos(bearing);
      double py = range * sin(bearing);
      x_ << px, py, 0.0, 0.0, 0.0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;
    }
    P_.fill(0.0);
    P_(0, 0) = .04;//.03
    P_(1, 1) = .06;
    P_(2, 2) = 10;//0.6
    P_(3, 3) = 10;//0.02
    P_(4, 4) = 10;//0.02
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  } else {
    /*****************************************************************************
       *  Prediction
       ****************************************************************************/

      //compute the time elapsed between the current and previous measurements
      double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
      if (dt < NEAR_ZERO) {
        return;
      }
      time_us_ = meas_package.timestamp_;

      Prediction(dt);

      /*****************************************************************************
       *  Update
       ****************************************************************************/

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
      } else {
        UpdateLidar(meas_package);
      }

      // print the output
      cout << "x_ = " << x_ << endl;
      cout << "P_ = " << P_ << endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd Xsig_aug1 = AugmentedSigmaPoints(Xsig_aug);
  SigmaPointPrediction(Xsig_aug1, delta_t);
  PredictMeanAndCovariance();
}

MatrixXd UKF::AugmentedSigmaPoints(MatrixXd &Xsig_aug) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  return Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t) {

    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
      //extract values for better readability
      double p_x = Xsig_aug(0,i);
      double p_y = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double yaw = Xsig_aug(3,i);
      double yawd = Xsig_aug(4,i);
      double nu_a = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);

      //predicted state values
      double px_p, py_p;

      //avoid division by zero
      if (fabs(yawd) > NEAR_ZERO) {
          px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
          py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
      }
      else {
          px_p = p_x + v*delta_t*cos(yaw);
          py_p = p_y + v*delta_t*sin(yaw);
      }

      double v_p = v;
      double yaw_p = yaw + yawd*delta_t;
      double yawd_p = yawd;

      //add noise
      px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
      py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
      v_p = v_p + nu_a*delta_t;

      yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
      yawd_p = yawd_p + nu_yawdd*delta_t;

      //write predicted sigma point into right column
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
    }
}

void UKF::PredictMeanAndCovariance() {
  //predicted state mean
  //x_.fill(0.0);
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = ConvertToPiRange(x_diff(3));


    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  cout << "Predicted state" << endl;
  cout << x_ << endl;
  cout << "Predicted covariance matrix" << endl;
  cout << P_ << endl;
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
  int n_z = 2;
  MatrixXd H = MatrixXd(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;
  MatrixXd Ht = H.transpose();

  VectorXd z_pred = H * x_;
  VectorXd y = meas_package.raw_measurements_ - z_pred;

  MatrixXd S = H * P_ * Ht;

  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  S = S + R;

  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  x_ = x_ + (K * y);
  P_ = (I - K * H) * P_;
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
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  //innovation covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  PredictRadarMeasurement(Zsig, z_pred, S, n_z);
  UpdateRadarMeasurement(Zsig, z_pred, S, meas_package.raw_measurements_, n_z);
}

void UKF::PredictRadarMeasurement(MatrixXd &Zsig, VectorXd  &z_pred, MatrixXd &S, int n_z) {
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);

    if (p_x < NEAR_ZERO && p_y < NEAR_ZERO) {
      Zsig(1,i) = 0.0;
    } else {
      Zsig(1,i) = atan2(p_y,p_x);
    }

    if (Zsig(0, i) < NEAR_ZERO) {
      Zsig(2, i) = 0.0;
    } else {
      Zsig(2, i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);
    }
  }

cout << "Zsig: " << endl << Zsig << endl;

  //mean predicted measurement
  z_pred.fill(0.0);
  z_pred = Zsig * weights_;

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = ConvertToPiRange(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //print result
  cout << "Radar z_pred: " << endl << z_pred << endl;
  cout << "Radar S: " << endl << S << endl;
}

void UKF::UpdateRadarMeasurement(MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, VectorXd &z, int n_z) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = ConvertToPiRange(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = ConvertToPiRange(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = ConvertToPiRange(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  cout << "Updated state x: " << endl << x_ << endl;
  cout << "Updated state covariance P: " << endl << P_ << endl;

}
