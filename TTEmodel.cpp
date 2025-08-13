$PLUGIN mrgx Rcpp

// mod401scm_PFS_1224data
// mod402_AEDC_1224data
// mod402_G2GI_1224data

$PARAM

TVCL = 20.4
TVQ = 18.5
TVV1 = 4010
TVV2 = 2450
BASE = exp(-9.19)
SHP = 1.52
EC50_PFS = 0.03
EMAX_PFS = 3
CAVGC1 = -1
CL = -1
V1 = -1
Q = -1
V2 = -1
WT = -1

$SET delta=24, end=6000

$CMT CENT PERI AUC HAZPFS HAZ


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
double DEL = 1E-6;
double DRUGEF_PFS = pow(CAVGC1/2.6,0.9);

$ODE

if(IPRED > CMAX) maxt = SOLVERTIME;
if(IPRED > CMAX) CMAX = IPRED;
if(IPRED < CMIN & SOLVERTIME > 1) CMIN = IPRED;
dxdt_CENT = -CL*CENT/V1 - Q*CENT/V1 + Q*PERI/V2;
dxdt_PERI = -(-Q*CENT/V1 + Q*PERI/V2);
dxdt_AUC = CENT/V1;
double IPRED = CENT/V1;
double CAVG = AUC/(SOLVERTIME+DEL);
dxdt_HAZPFS = (BASE*SHP*pow((BASE*(SOLVERTIME+DEL)),(SHP-1)))*(DRUGEF_PFS);
dxdt_HAZ = (BASE*SHP*pow((BASE*(SOLVERTIME+DEL)),(SHP-1)));

$TABLE
capture CP = CENT/V1;
capture CMHZ_PFS  = HAZPFS ;
capture ST_PFS    = exp(-CMHZ_PFS)   ;
capture CMHZ  = HAZ ;
capture ST    = exp(-CMHZ)   ;

$CAPTURE IPRED CMAX CMIN CAVG
DRUGEF_PFS WT
CL V1 Q V2
