#ifndef _ANA_CONSTANTS_H
#define _ANA_CONSTANTS_H	
//global variables and binning

string posCut = "atan(y1_st3/y1_st1)>1.086 && atan(y1_st3/y1_st1) <1.562&& z1_v > -45.0  && z1_v < 125.0&& x1_t*x1_t+(y1_t-1.6)*(y1_t-1.6) <1425.0 && x1_d*x1_d+(y1_d -1.6)*(y1_d-1.6) <141.0&& abs(px1_st1-px1_st3)> 0.413 &&  abs(px1_st1-px1_st3)<0.420 && abs(py1_st1-py1_st3)< 0.005&& abs(pz1_st1-pz1_st3)<0.074&&(chisq1_target-chisq1_dump) >-2.5 && (chisq1_target-chisq1_dump)<77.5&& pz1_st1 > 10 && pz1_st1 < 78 && nHits1 > 13 &&(y1_st1*y1_st3) > 0 && chisq1_dump/chisq1_upstream < 0.165&& nHits1St1 >4 &&";

string negCut ="atan(y2_st3/y2_st1)>0.666 && atan(y2_st3/y2_st1) <1.506&& z2_v > -45.0 && z2_v < 135.0&& x2_t*x2_t+(y2_t-1.6)*(y2_t-1.6) <1485.0 && x2_d*x2_d+(y2_d-1.6)*(y2_d-1.6) <201.0&& abs(px2_st1-px2_st3)> 0.414 &&  abs(px2_st1-px2_st3)<0.420 && abs(py2_st1-py2_st3)< 0.005&& abs(pz2_st1-pz2_st3)<0.05&&(chisq2_target-chisq2_dump) >-2.5 && (chisq2_target-chisq2_dump)<77.5&&  pz2_st1 > 9 && pz2_st1 < 78 && nHits2 > 13 &&(y2_st1*y2_st3) > 0&& chisq2_dump/chisq2_upstream < 0.195 && nHits2St1 >4 &&";

string dimuCut= "dz>0.0 && dz< 125.0 && abs(dx) < 0.262 && abs(dy-1.6) < 0.225 &&  abs(dpx) < 1.675 && abs(dpy) < 2.4 && dpz > 40.0 && dpz < 114.0 && abs(x1_st1 + x2_st1) < 30 && dpx*dpx + dpy*dpy < 5.3 && dx*dx + (dy-1.6)*(dy-1.6) < 0.105  && y1_st3*y2_st3 < 0 && nHits1 + nHits2 > 30 && nHits1St1 + nHits2St1 >9 &&chisq_dimuon>0.0&&chisq_dimuon<11.5&&(chisq1_dump+chisq2_dump-chisq_dimuon)<19.27&&abs(z1_v - z2_v)<125.0&& abs(trackSeparation) < 270 &&";
string physics = "xF > -0.18 && xF <0.94&& xT > 0.007 && xT <0.54 && mass>5.0 && mass <8.0 && pT < 2.0 && D1<400";
string trigger = "fpga1==1&&";
//============Final Cut
string cutTrue ="true_mass>=4.0";
string cutRecoData =posCut+negCut+dimuCut+physics+ "&&D1<400&& D2<400 && D3 <400"; 
string cutRecoMC= cutTrue+ "&&"+trigger +cutRecoData;
#endif
