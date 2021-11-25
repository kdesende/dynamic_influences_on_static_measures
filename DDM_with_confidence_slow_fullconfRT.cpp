// in r: install rcpp library and ZigguratRCPP library, then source this script and the function becomes available (source('DDM_KSfitting'))
#include <Rcpp.h>
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;
// [[Rcpp::export]]
NumericMatrix DDM_with_confidence_slow_fullconfRT(double v, double a, double ter, double z, int ntrials, double s, double dt, NumericVector t2distribution, double postdriftmod) {

  // initialize output
  NumericMatrix DATA(ntrials,6);

  // loop over trials
  for (int i = 0; i < ntrials; i++) {

    // initalize variables
    int resp = -1;
    int acc = 0;
    double t2time = 0;
    double evidence = a*z;
    double t = 0;

    // Decisional processing
    while (t > -1){

      t = t + dt;
      evidence = evidence + v * dt + s * sqrt(dt) * zigg.norm();

      if (evidence >= a){
        resp = 1;
        evidence=a;
        if (v > 0){
          acc = 1;
        }
        break;
      } else if (evidence <= 0) {
        resp = -1;
        evidence=0;
        if (v < 0){
          acc = 1;
        }
        break;
      }
    }

    DATA(i,0) = (t + ter);
    DATA(i,1) = resp;
    DATA(i,2) = acc;
    
    t2time = t2distribution(i);

    //Post-decisional processing
    double v_post = v * postdriftmod;
    for (int j = 0; j < t2time/dt; j++){
      t = t + dt;
      evidence = evidence + v_post * dt + s * sqrt(dt) * zigg.norm();
    }
    DATA(i,3) = evidence;
    DATA(i,4) = (t + ter);
    if (resp == 1){
      DATA(i,5) = evidence-a;
    }
    if(resp == -1){
      DATA(i,5) = -evidence;
    }
  }

  return DATA; //RT, resp,accuracy, evidence2, rt2, confidence
}
