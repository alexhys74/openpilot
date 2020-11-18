/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3097991058676388603);
void inv_err_fun(double *nom_x, double *true_x, double *out_6781449444750105927);
void H_mod_fun(double *state, double *out_8275655445109615529);
void f_fun(double *state, double dt, double *out_4520924155593214072);
void F_fun(double *state, double dt, double *out_7526409145971270572);
void h_3(double *state, double *unused, double *out_4400580297741421371);
void H_3(double *state, double *unused, double *out_4263837380245414145);
void h_4(double *state, double *unused, double *out_7510251667274556838);
void H_4(double *state, double *unused, double *out_3462043278788988877);
void h_9(double *state, double *unused, double *out_3261028494164123787);
void H_9(double *state, double *unused, double *out_5828164258346798549);
void h_10(double *state, double *unused, double *out_672471820299075556);
void H_10(double *state, double *unused, double *out_6911236767292881296);
void h_12(double *state, double *unused, double *out_6711799318218590937);
void H_12(double *state, double *unused, double *out_2033519118481025477);
void h_13(double *state, double *unused, double *out_4336546972390840157);
void H_13(double *state, double *unused, double *out_8758755579235800108);
void h_14(double *state, double *unused, double *out_3261028494164123787);
void H_14(double *state, double *unused, double *out_5828164258346798549);
void h_19(double *state, double *unused, double *out_3907083411332888491);
void H_19(double *state, double *unused, double *out_5780489044463970239);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);