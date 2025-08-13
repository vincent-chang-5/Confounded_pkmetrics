$PLUGIN mrgx Rcpp

// mod401scm_PFS_1224data
// mod402_AEDC_1224data
// mod402_G2GI_1224data

$PARAM

TVCL = 20.4
TVQ = 18.5
TVV1 = 4010
TVV2 = 2450
BASE = exp(-8.19)
SHP = 1.52
EC50_PFS = 0.03
EMAX_PFS = 3

$SET delta=24, end=6000

$CMT CENT PERI AUC


$GLOBAL
double maxt = 0;
double CMAX = 0;
double CMIN = 100000;


$MAIN
if(NEWIND<=1) {
  maxt = 0;
  CMAX = 0;
  CMIN = 100000000;
}
if(NEWIND <=1) {
  double WT = mrgx::rnorm(85,19.3,40,150);
}
double CL = TVCL*exp(ECL)*pow((WT/70), 0.75);
double V1 = TVV1*exp(EV1)*pow((WT/70), 1);
double Q = TVQ*pow((WT/70), 0.75);
double V2 = TVV2*pow((WT/70), 1);
double DEL = 1E-6;

$OMEGA @name IIV @block @labels ECL EV1
0.0734 0.0525 0.133

$ODE

if(IPRED > CMAX) maxt = SOLVERTIME;
if(IPRED > CMAX) CMAX = IPRED;
if(IPRED < CMIN & SOLVERTIME > 1) CMIN = IPRED;
dxdt_CENT = -CL*CENT/V1 - Q*CENT/V1 + Q*PERI/V2;
dxdt_PERI = -(-Q*CENT/V1 + Q*PERI/V2);
dxdt_AUC = CENT/V1;
double IPRED = CENT/V1;
double CAVG = AUC/(SOLVERTIME+DEL);

$CAPTURE maxt IPRED CMAX CMIN ECL EV1 CAVG
WT
CL V1 Q V2
