
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8007853367845348407) {
   out_8007853367845348407[0] = delta_x[0] + nom_x[0];
   out_8007853367845348407[1] = delta_x[1] + nom_x[1];
   out_8007853367845348407[2] = delta_x[2] + nom_x[2];
   out_8007853367845348407[3] = delta_x[3] + nom_x[3];
   out_8007853367845348407[4] = delta_x[4] + nom_x[4];
   out_8007853367845348407[5] = delta_x[5] + nom_x[5];
   out_8007853367845348407[6] = delta_x[6] + nom_x[6];
   out_8007853367845348407[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5980494832453821423) {
   out_5980494832453821423[0] = -nom_x[0] + true_x[0];
   out_5980494832453821423[1] = -nom_x[1] + true_x[1];
   out_5980494832453821423[2] = -nom_x[2] + true_x[2];
   out_5980494832453821423[3] = -nom_x[3] + true_x[3];
   out_5980494832453821423[4] = -nom_x[4] + true_x[4];
   out_5980494832453821423[5] = -nom_x[5] + true_x[5];
   out_5980494832453821423[6] = -nom_x[6] + true_x[6];
   out_5980494832453821423[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5970333948667268010) {
   out_5970333948667268010[0] = 1.0;
   out_5970333948667268010[1] = 0.0;
   out_5970333948667268010[2] = 0.0;
   out_5970333948667268010[3] = 0.0;
   out_5970333948667268010[4] = 0.0;
   out_5970333948667268010[5] = 0.0;
   out_5970333948667268010[6] = 0.0;
   out_5970333948667268010[7] = 0.0;
   out_5970333948667268010[8] = 0.0;
   out_5970333948667268010[9] = 1.0;
   out_5970333948667268010[10] = 0.0;
   out_5970333948667268010[11] = 0.0;
   out_5970333948667268010[12] = 0.0;
   out_5970333948667268010[13] = 0.0;
   out_5970333948667268010[14] = 0.0;
   out_5970333948667268010[15] = 0.0;
   out_5970333948667268010[16] = 0.0;
   out_5970333948667268010[17] = 0.0;
   out_5970333948667268010[18] = 1.0;
   out_5970333948667268010[19] = 0.0;
   out_5970333948667268010[20] = 0.0;
   out_5970333948667268010[21] = 0.0;
   out_5970333948667268010[22] = 0.0;
   out_5970333948667268010[23] = 0.0;
   out_5970333948667268010[24] = 0.0;
   out_5970333948667268010[25] = 0.0;
   out_5970333948667268010[26] = 0.0;
   out_5970333948667268010[27] = 1.0;
   out_5970333948667268010[28] = 0.0;
   out_5970333948667268010[29] = 0.0;
   out_5970333948667268010[30] = 0.0;
   out_5970333948667268010[31] = 0.0;
   out_5970333948667268010[32] = 0.0;
   out_5970333948667268010[33] = 0.0;
   out_5970333948667268010[34] = 0.0;
   out_5970333948667268010[35] = 0.0;
   out_5970333948667268010[36] = 1.0;
   out_5970333948667268010[37] = 0.0;
   out_5970333948667268010[38] = 0.0;
   out_5970333948667268010[39] = 0.0;
   out_5970333948667268010[40] = 0.0;
   out_5970333948667268010[41] = 0.0;
   out_5970333948667268010[42] = 0.0;
   out_5970333948667268010[43] = 0.0;
   out_5970333948667268010[44] = 0.0;
   out_5970333948667268010[45] = 1.0;
   out_5970333948667268010[46] = 0.0;
   out_5970333948667268010[47] = 0.0;
   out_5970333948667268010[48] = 0.0;
   out_5970333948667268010[49] = 0.0;
   out_5970333948667268010[50] = 0.0;
   out_5970333948667268010[51] = 0.0;
   out_5970333948667268010[52] = 0.0;
   out_5970333948667268010[53] = 0.0;
   out_5970333948667268010[54] = 1.0;
   out_5970333948667268010[55] = 0.0;
   out_5970333948667268010[56] = 0.0;
   out_5970333948667268010[57] = 0.0;
   out_5970333948667268010[58] = 0.0;
   out_5970333948667268010[59] = 0.0;
   out_5970333948667268010[60] = 0.0;
   out_5970333948667268010[61] = 0.0;
   out_5970333948667268010[62] = 0.0;
   out_5970333948667268010[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6337864808026530212) {
   out_6337864808026530212[0] = state[0];
   out_6337864808026530212[1] = state[1];
   out_6337864808026530212[2] = state[2];
   out_6337864808026530212[3] = state[3];
   out_6337864808026530212[4] = state[4];
   out_6337864808026530212[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6337864808026530212[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6337864808026530212[7] = state[7];
}
void F_fun(double *state, double dt, double *out_5220790824763354224) {
   out_5220790824763354224[0] = 1;
   out_5220790824763354224[1] = 0;
   out_5220790824763354224[2] = 0;
   out_5220790824763354224[3] = 0;
   out_5220790824763354224[4] = 0;
   out_5220790824763354224[5] = 0;
   out_5220790824763354224[6] = 0;
   out_5220790824763354224[7] = 0;
   out_5220790824763354224[8] = 0;
   out_5220790824763354224[9] = 1;
   out_5220790824763354224[10] = 0;
   out_5220790824763354224[11] = 0;
   out_5220790824763354224[12] = 0;
   out_5220790824763354224[13] = 0;
   out_5220790824763354224[14] = 0;
   out_5220790824763354224[15] = 0;
   out_5220790824763354224[16] = 0;
   out_5220790824763354224[17] = 0;
   out_5220790824763354224[18] = 1;
   out_5220790824763354224[19] = 0;
   out_5220790824763354224[20] = 0;
   out_5220790824763354224[21] = 0;
   out_5220790824763354224[22] = 0;
   out_5220790824763354224[23] = 0;
   out_5220790824763354224[24] = 0;
   out_5220790824763354224[25] = 0;
   out_5220790824763354224[26] = 0;
   out_5220790824763354224[27] = 1;
   out_5220790824763354224[28] = 0;
   out_5220790824763354224[29] = 0;
   out_5220790824763354224[30] = 0;
   out_5220790824763354224[31] = 0;
   out_5220790824763354224[32] = 0;
   out_5220790824763354224[33] = 0;
   out_5220790824763354224[34] = 0;
   out_5220790824763354224[35] = 0;
   out_5220790824763354224[36] = 1;
   out_5220790824763354224[37] = 0;
   out_5220790824763354224[38] = 0;
   out_5220790824763354224[39] = 0;
   out_5220790824763354224[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5220790824763354224[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5220790824763354224[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5220790824763354224[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5220790824763354224[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5220790824763354224[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5220790824763354224[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5220790824763354224[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5220790824763354224[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5220790824763354224[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5220790824763354224[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5220790824763354224[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5220790824763354224[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5220790824763354224[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5220790824763354224[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5220790824763354224[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5220790824763354224[56] = 0;
   out_5220790824763354224[57] = 0;
   out_5220790824763354224[58] = 0;
   out_5220790824763354224[59] = 0;
   out_5220790824763354224[60] = 0;
   out_5220790824763354224[61] = 0;
   out_5220790824763354224[62] = 0;
   out_5220790824763354224[63] = 1;
}
void h_25(double *state, double *unused, double *out_1295897918763796329) {
   out_1295897918763796329[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8617664870809116450) {
   out_8617664870809116450[0] = 0;
   out_8617664870809116450[1] = 0;
   out_8617664870809116450[2] = 0;
   out_8617664870809116450[3] = 0;
   out_8617664870809116450[4] = 0;
   out_8617664870809116450[5] = 0;
   out_8617664870809116450[6] = 1;
   out_8617664870809116450[7] = 0;
}
void h_24(double *state, double *unused, double *out_5308920038595990296) {
   out_5308920038595990296[0] = state[4];
   out_5308920038595990296[1] = state[5];
}
void H_24(double *state, double *unused, double *out_610989344707219134) {
   out_610989344707219134[0] = 0;
   out_610989344707219134[1] = 0;
   out_610989344707219134[2] = 0;
   out_610989344707219134[3] = 0;
   out_610989344707219134[4] = 1;
   out_610989344707219134[5] = 0;
   out_610989344707219134[6] = 0;
   out_610989344707219134[7] = 0;
   out_610989344707219134[8] = 0;
   out_610989344707219134[9] = 0;
   out_610989344707219134[10] = 0;
   out_610989344707219134[11] = 0;
   out_610989344707219134[12] = 0;
   out_610989344707219134[13] = 1;
   out_610989344707219134[14] = 0;
   out_610989344707219134[15] = 0;
}
void h_30(double *state, double *unused, double *out_9214349337968772811) {
   out_9214349337968772811[0] = state[4];
}
void H_30(double *state, double *unused, double *out_545790703131247336) {
   out_545790703131247336[0] = 0;
   out_545790703131247336[1] = 0;
   out_545790703131247336[2] = 0;
   out_545790703131247336[3] = 0;
   out_545790703131247336[4] = 1;
   out_545790703131247336[5] = 0;
   out_545790703131247336[6] = 0;
   out_545790703131247336[7] = 0;
}
void h_26(double *state, double *unused, double *out_7619273719607758096) {
   out_7619273719607758096[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8672107487314229750) {
   out_8672107487314229750[0] = 0;
   out_8672107487314229750[1] = 0;
   out_8672107487314229750[2] = 0;
   out_8672107487314229750[3] = 0;
   out_8672107487314229750[4] = 0;
   out_8672107487314229750[5] = 0;
   out_8672107487314229750[6] = 0;
   out_8672107487314229750[7] = 1;
}
void h_27(double *state, double *unused, double *out_9137457230113761564) {
   out_9137457230113761564[0] = state[3];
}
void H_27(double *state, double *unused, double *out_788690960601286298) {
   out_788690960601286298[0] = 0;
   out_788690960601286298[1] = 0;
   out_788690960601286298[2] = 0;
   out_788690960601286298[3] = 1;
   out_788690960601286298[4] = 0;
   out_788690960601286298[5] = 0;
   out_788690960601286298[6] = 0;
   out_788690960601286298[7] = 0;
}
void h_29(double *state, double *unused, double *out_1869446293172429638) {
   out_1869446293172429638[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3394270571481645070) {
   out_3394270571481645070[0] = 0;
   out_3394270571481645070[1] = 1;
   out_3394270571481645070[2] = 0;
   out_3394270571481645070[3] = 0;
   out_3394270571481645070[4] = 0;
   out_3394270571481645070[5] = 0;
   out_3394270571481645070[6] = 0;
   out_3394270571481645070[7] = 0;
}
void h_28(double *state, double *unused, double *out_3849094725402729368) {
   out_3849094725402729368[0] = state[5];
   out_3849094725402729368[1] = state[6];
}
void H_28(double *state, double *unused, double *out_412651612576102076) {
   out_412651612576102076[0] = 0;
   out_412651612576102076[1] = 0;
   out_412651612576102076[2] = 0;
   out_412651612576102076[3] = 0;
   out_412651612576102076[4] = 0;
   out_412651612576102076[5] = 1;
   out_412651612576102076[6] = 0;
   out_412651612576102076[7] = 0;
   out_412651612576102076[8] = 0;
   out_412651612576102076[9] = 0;
   out_412651612576102076[10] = 0;
   out_412651612576102076[11] = 0;
   out_412651612576102076[12] = 0;
   out_412651612576102076[13] = 0;
   out_412651612576102076[14] = 1;
   out_412651612576102076[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
