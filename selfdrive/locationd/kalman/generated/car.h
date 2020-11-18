/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8007853367845348407);
void inv_err_fun(double *nom_x, double *true_x, double *out_5980494832453821423);
void H_mod_fun(double *state, double *out_5970333948667268010);
void f_fun(double *state, double dt, double *out_6337864808026530212);
void F_fun(double *state, double dt, double *out_5220790824763354224);
void h_25(double *state, double *unused, double *out_1295897918763796329);
void H_25(double *state, double *unused, double *out_8617664870809116450);
void h_24(double *state, double *unused, double *out_5308920038595990296);
void H_24(double *state, double *unused, double *out_610989344707219134);
void h_30(double *state, double *unused, double *out_9214349337968772811);
void H_30(double *state, double *unused, double *out_545790703131247336);
void h_26(double *state, double *unused, double *out_7619273719607758096);
void H_26(double *state, double *unused, double *out_8672107487314229750);
void h_27(double *state, double *unused, double *out_9137457230113761564);
void H_27(double *state, double *unused, double *out_788690960601286298);
void h_29(double *state, double *unused, double *out_1869446293172429638);
void H_29(double *state, double *unused, double *out_3394270571481645070);
void h_28(double *state, double *unused, double *out_3849094725402729368);
void H_28(double *state, double *unused, double *out_412651612576102076);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
