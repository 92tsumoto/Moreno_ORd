#include "syspara.h"

void comp_ina(double x[])
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	//printf("iV=%d,V1=%f,V2=%f,d1=%f,d2=%f\n",iV,V1,V2,d1,d2);

	// Fast Na current

	/* O'Hara-Rudy Original INa current
	ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
	ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
	ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
	ina.tauh_fast = ina.Ttauh_fast[iV]*d2 + ina.Ttauh_fast[iV+1]*d1;
	ina.tauh_slow = ina.Ttauh_slow[iV]*d2 + ina.Ttauh_slow[iV+1]*d1;
	ina.h = ina.Ah_fast*x[2] + ina.Ah_slow*x[3];
	ina.jss = ina.hss;
	ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;

	ina.jss = ina.hss;
	ina.tauj = 2.038+1.0/(0.02136*exp(-(x[0]+100.6)/8.281)+0.3052*exp((x[0]+0.9941)/38.45));
	
	ina.hCaMKss = ina.ThCaMKss[iV]*d2 + ina.ThCaMKss[iV+1]*d1;
	ina.tauh_CaMK_slow = 3.0*ina.tauh_slow;
	ina.hCaMK = ina.Ah_CaMK_fast*x[2] + ina.Ah_CaMK_slow*x[5];	// h_CaMK_fast = h_fast (x[2])
	ina.hCaMK = ina.Ah_fast*x[2] + ina.Ah_slow*x[5];	// h_CaMK_fast = h_fast (x[2])

	ina.jCaMKss = ina.jss; //jss=hss;
	ina.tauj_CaMK = 1.46*ina.tauj;
	
	ina.fast_pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);

	ina.fast = ina.Gna_fast*(x[0]-var.Ena)*x[1]*x[1]*x[1]*((1.0-ina.fast_pCaMK)*ina.h*x[4]+ina.fast_pCaMK*ina.hCaMK*x[6]);
	*/
	
// Fast Sodium Current (time dependant) */

	ina.kd_open = ina.kd0*exp( (-0.7*x[0]*F)/(R*T));
	ina.koff = ina.kd_open * ina.diffusion;
	ina.kcoff = ina.koff;
	// Drug free part
		ina.a11 = ina.Tfac*(ina.Ta11[iV]*d2 + ina.Ta11[iV+1]*d1);
		ina.a12 = ina.Tfac*(ina.Ta12[iV]*d2 + ina.Ta12[iV+1]*d1);
		ina.a13 = ina.Tfac*(ina.Ta13[iV]*d2 + ina.Ta13[iV+1]*d1);
		ina.b11 = ina.Tfac*(ina.Tb11[iV]*d2 + ina.Tb11[iV+1]*d1);
		ina.b12 = ina.Tfac*(ina.Tb12[iV]*d2 + ina.Tb12[iV+1]*d1);
		ina.b13 = ina.Tfac*(ina.Tb13[iV]*d2 + ina.Tb13[iV+1]*d1);
		ina.a3  = ina.Tfac*(ina.Ta3[iV]*d2 + ina.Ta3[iV+1]*d1);
		ina.b3  = ina.Tfac*(ina.Tb3[iV]*d2 + ina.Tb3[iV+1]*d1);
		ina.a2  = ina.Tfac*(ina.Ta2[iV]*d2 + ina.Ta2[iV+1]*d1);
		ina.b2  = (ina.a13*ina.a2*ina.a3)/(ina.b13*ina.b3);
		ina.ax  = 3.4229E-2*ina.a2;
		ina.bx  = 1.7898E-2*ina.a3;

	if(var.drug_type==0){	// Flecanide
		ina.ax1 = 5.7839E-5*ina.ax;
		ina.bx1 = 1.6689E-8*ina.bx;
		ina.a13c = 3.6324E-3*ina.a13;
		ina.a22 = 1.4847E+3*ina.a2;
		ina.b33 = 1.7352E-6*ina.b3;
		ina.a33 = 6.7505E-5*ina.a3;
		ina.a44 = 2.4135*ina.a2;
		ina.b44 = 4.9001E-2*ina.a3;

		ina.ax2 = 2.6126E-1*ina.ax;
		ina.a13n = 2.6452*ina.a13;
		ina.a_22 = 4.2385E+1*ina.a2;
		ina.b_33 = 2.1181*ina.b3;
		ina.a_44 = 1.0326E-3*ina.a2;
		//ina.b_44 = 2.1378E-2*ina.a3;
		ina.b_44 = 2.1378E-2*ina.a2;

	} else if(var.drug_type == 1){	// Lidocaine
		ina.ax1 = 6.3992E-7*ina.ax;
		ina.bx1 = 1.3511*ina.bx;
		ina.a13c = 5.6974E-3*ina.a13;
		ina.a22 = 6.7067E-6*ina.a2;
		ina.b33 = 1.9698E-5*ina.b3;
		ina.a33 = 3.2976*ina.a3;
		ina.a44 = 0.0;
		ina.b44 = 0.0;

		ina.ax2 = 1.3110E-1*ina.ax;
		ina.a13n = 8.4559E+1*ina.a13;
		ina.a_22 = 1.7084E-5*ina.a2;
		ina.b_33 = 4.8477*ina.b3;
		ina.a_44 = 0.0;
		ina.b_44 = 0.0;
	}

		if(var.drug == 0.0 || ina.drug_charged == 0.0 ){
			ina.b13c = 0.0;
		} else {
			ina.b13c = (ina.b13*ina.kcon*ina.koff*ina.a13c)/(ina.kon*ina.kcoff*ina.a13);
		}
		if(var.drug == 0.0 || ina.drug_neutral == 0.0){
			ina.b13n = 0.0;
		} else {
			ina.b13n = (ina.b13*ina.kc_on*ina.a13n*ina.k_off)/(ina.kc_off*ina.a13*ina.k_on);
		}
		if(var.drug == 0.0 || ina.drug_neutral == 0.0){
			ina.bx2 = 0.0;
		} else {
			ina.bx2 = (ina.bx*ina.k_on*ina.ax2*ina.ki_off)/(ina.ax*ina.ki_on*ina.k_off);
		}
		if(ina.b13c == 0.0){
			ina.b22 = 0.0;
		} else {
			ina.b22 = (ina.a13c*ina.a22*ina.a33)/(ina.b13c*ina.b33);
		}
		if(ina.b13n == 0.0){
			ina.b_22 = 0.0;
		} else {
			ina.b_22 = (ina.a_33*ina.a13n*ina.a_22)/(ina.b_33*ina.b13n);
		}
		if(var.drug==0.0 || ina.drug_neutral == 0.0){
			ina.a_33 = 0.0;
		} else {
			ina.a_33 = (ina.ki_off*ina.a3*ina.kc_on*ina.b_33)/(ina.ki_on*ina.kc_off*ina.b3);
		}

	ina.fast = ina.Gna_fast*x[1]*(x[0] - var.Ena);

	// late Na current
	ina.mlss = ina.Tmlss[iV]*d2 + ina.Tmlss[iV+1]*d1;
	//ina.tauml = ina.taum;
	ina.tauml = ina.Ttauml[iV]*d2 + ina.Ttauml[iV+1]*d1;
	ina.hlss = ina.Thlss[iV]*d2 + ina.Thlss[iV+1]*d1;
	ina.hlCaMKss = ina.ThlCaMKss[iV]*d2 + ina.ThlCaMKss[iV+1]*d1;

	ina.late_pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);
	ina.late = ina.Gna_late*(x[0]-var.Ena)*x[27]*((1.0-ina.late_pCaMK)*x[28] + ina.late_pCaMK*x[29]);
	
	ina.total = ina.fast + ina.late;
}

// Ito Transient Outward Current
void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.ass = ito.Tass[iV]*d2 + ito.Tass[iV+1]*d1;
	ito.taua = ito.Ttaua[iV]*d2 + ito.Ttaua[iV+1]*d1;

	ito.iss = ito.Tiss[iV]*d2 + ito.Tiss[iV+1]*d1;
	if(var.celltype==1){
		ito.depi = ito.Tdepi[iV]*d2 + ito.Tdepi[iV+1]*d1;
	} else {
		ito.depi = 1.0;
	}
	ito.taui_fast = ito.depi*(ito.Ttaui_fast[iV]*d2 + ito.Ttaui_fast[iV+1]*d1);
	ito.taui_slow = ito.depi*(ito.Ttaui_slow[iV]*d2 + ito.Ttaui_slow[iV+1]*d1);
	ito.Ai_fast = ito.TAi_fast[iV]*d2 + ito.TAi_fast[iV+1]*d1;
	ito.Ai_slow = 1.0 - ito.Ai_fast;
	ito.i = ito.Ai_fast*x[31] + ito.Ai_slow*x[32];

	ito.aCaMKss = ito.TaCaMKss[iV]*d2 + ito.TaCaMKss[iV+1]*d1;
	//ito.taua_CaMK = ito.taua;

	//ito.iCaMKss = ito.iss;
	ito.deltaCaMK_dev = ito.TdeltaCaMK_dev[iV]*d2 + ito.TdeltaCaMK_dev[iV+1]*d1;
	ito.deltaCaMK_rec = ito.TdeltaCaMK_rec[iV]*d2 + ito.TdeltaCaMK_rec[iV+1]*d1;
	ito.taui_CaMK_fast = ito.taui_fast*ito.deltaCaMK_dev*ito.deltaCaMK_rec;
	ito.taui_CaMK_slow = ito.taui_slow*ito.deltaCaMK_dev*ito.deltaCaMK_rec;
	//ito.Ai_CaMK_fast = ito.Ai_fast;
	//ito.Ai_CaMK_slow = ito.Ai_slow;
	//ito.iCaMK = ito.Ai_CaMK_fast*x[34] + ito.Ai_CaMK_slow*x[35];
	ito.iCaMK = ito.Ai_fast*x[34] + ito.Ai_slow*x[35];

	ito.pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);
	ito.ik = ito.Gto*(x[0]-var.Ek)*((1.0-ito.pCaMK)*x[30]*ito.i + ito.pCaMK*x[33]*ito.iCaMK);

}


// L-type calcium current
void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

// VDA
	ical.dss = ical.Tdss[iV]*d2 + ical.Tdss[iV+1]*d1;
	ical.taud = ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1;
// VDI 
	ical.fss = ical.Tfss[iV]*d2 + ical.Tfss[iV+1]*d1;
	ical.tauf_fast = ical.Ttauf_fast[iV]*d2 + ical.Ttauf_fast[iV+1]*d1;
	ical.tauf_slow = ical.Ttauf_slow[iV]*d2 + ical.Ttauf_slow[iV+1]*d1;
	ical.f = ical.Af_fast*x[37] + ical.Af_slow*x[38];

// CDI 
	//ical.fcass = ical.fss;
	ical.taufca_fast = ical.Ttaufca_fast[iV]*d2 + ical.Ttaufca_fast[iV+1]*d1;
	ical.taufca_slow = ical.Ttaufca_slow[iV]*d2 + ical.Ttaufca_slow[iV+1]*d1;
	ical.Afca_fast = ical.TAfca_fast[iV]*d2 + ical.TAfca_fast[iV+1]*d1;
	ical.Afca_slow = 1.0 - ical.Afca_fast;
	ical.fca = ical.Afca_fast*x[39] + ical.Afca_slow*x[40];

	//ical.jcass = ical.fcass;
	//ical.taujca = 75.0; // (ms).

// CaMK(VDI)
    //ical.f_CaMKss = ical.fss;
	ical.tauf_CaMK_fast = 2.5*ical.tauf_fast;
	//ical.Af_CaMK_fast = ical.Af_fast;
	//ical.Af_CaMK_slow = ical.Af_slow;
	//ical.f_CaMK_slow = x[38]; //(x[38] = f_slow )
	ical.f_CaMK = ical.Af_fast*x[42] + ical.Af_slow*x[38];

// CaMK(CDI)
	
	//ical.fca_CaMKss = ical.fss;
	ical.taufca_CaMK_fast = 2.5*ical.taufca_fast;
	//ical.Afca_CaMK_fast = ical.Afca_fast;
	//ical.Afca_CaMK_slow = ical.Afca_slow;
	//var.fca_CaMK_slow = x[40]; //(x[40] = fca_slow )
	ical.fca_CaMK = ical.Afca_fast*x[43] + ical.Afca_slow*x[40];

	//ical.kmn = 0.002;
	//ical.kp2n = 1000.0;
	ical.km2n = 1.0*x[41]; //(jca: recovery from CDI for ICaL )
	ical.tmp = (1.0+ical.kmn/x[58])*(1.0+ical.kmn/x[58])*(1.0+ical.kmn/x[58])*(1.0+ical.kmn/x[58]);
	ical.alpha_n = 1.0/((ical.kp2n/ical.km2n)+ical.tmp);
	
	ical.exp_Ca = ical.Texp_Ca[iV]*d2 + ical.Texp_Ca[iV+1]*d1;
	ical.exp_Na = ical.Texp_Na[iV]*d2 + ical.Texp_Na[iV+1]*d1;
	ical.exp_K = ical.Texp_K[iV]*d2 + ical.Texp_K[iV+1]*d1;

	//if(fabs(x[0])>1E-8){
		ical.phi_ca = (4.0*F*F*x[0]/R/T)*(ical.gacai*x[58]*ical.exp_Ca-ical.gacao*var.cao)/(ical.exp_Ca-1.0);
		ical.phi_na = (1.0*F*F*x[0]/R/T)*(ical.ganai*x[54]*ical.exp_Na-ical.ganao*var.nao)/(ical.exp_Na-1.0);
		ical.phi_k = (1.0*F*F*x[0]/R/T)*(ical.gaki*x[56]*ical.exp_K-ical.gako*var.ko)/(ical.exp_K-1.0);
	//} else {
	//	ical.phi_ca = 0.0;//(4.0*F*F*x[0]/R/T)*(ical.gacai*x[58]*ical.exp_Ca-ical.gacao*var.cao)/(ical.exp_Ca-1.0);
	//	ical.phi_na = 0.0;//(1.0*F*F*x[0]/R/T)*(ical.ganai*x[54]*ical.exp_Na-ical.ganao*var.nao)/(ical.exp_Na-1.0);
	//	ical.phi_k = 0.0;//(1.0*F*F*x[0]/R/T)*(ical.gaki*x[56]*ical.exp_K-ical.gako*var.ko)/(ical.exp_K-1.0);
	//}
	ical.ibarcal = ical.pca*ical.phi_ca;
	ical.ibarcana = ical.pcana*ical.phi_na;
	ical.ibarcak = ical.pcak*ical.phi_k;

	ical.ibarcal_CaMK = ical.pca_CaMK*ical.phi_ca;
	ical.ibarcana_CaMK = ical.pcana_CaMK*ical.phi_na;
	ical.ibarcak_CaMK = ical.pcak_CaMK*ical.phi_k;

	ical.pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);

	ical.ica =ical.ibarcal*x[36]*(1.0-ical.pCaMK)*(ical.f*(1.0-x[44])+ical.fca*x[44]*x[41])
				+ical.ibarcal_CaMK*x[36]*ical.pCaMK*(ical.f_CaMK*(1.0-x[44])+ical.fca_CaMK*x[44]*x[41]);
	ical.icana = ical.ibarcana*x[36]*(1.0- ical.pCaMK)*(ical.f*(1.0-x[44])+ical.fca*x[44]*x[41])
				+ ical.ibarcana_CaMK*x[36]*ical.pCaMK*(ical.f_CaMK*(1.0-x[44])+ical.fca_CaMK*x[44]*x[41]);
	ical.icak = ical.ibarcak*x[36]*(1.0- ical.pCaMK)*(ical.f*(1.0-x[44])+ical.fca*x[44]*x[41])
				+ ical.ibarcak_CaMK*x[36]*ical.pCaMK*(ical.f_CaMK*(1.0-x[44])+ical.fca_CaMK*x[44]*x[41]);

}

// Rapidly Activating Potassium Current 
void comp_ikr (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikr.xrss = ikr.Txrss[iV]*d2 + ikr.Txrss[iV+1]*d1;
	ikr.tauxr_fast = ikr.Ttauxr_fast[iV]*d2 + ikr.Ttauxr_fast[iV+1]*d1;
	ikr.tauxr_slow = ikr.Ttauxr_slow[iV]*d2 + ikr.Ttauxr_slow[iV+1]*d1;
	ikr.Axr_fast = ikr.TAxr_fast[iV]*d2 + ikr.TAxr_fast[iV+1]*d1;
	ikr.Axr_slow = 1.0 - ikr.Axr_fast;
	ikr.rkr = ikr.Trkr[iV]*d2 + ikr.Trkr[iV+1]*d1;

	ikr.xr = ikr.Axr_fast*x[45] + ikr.Axr_slow*x[46];

	ikr.ik = ikr.Gkr*ikr.rategkr*ikr.xr*ikr.rkr*(x[0]-var.Ek);
	//ikr.ik = 0.15*ikr.Gkr*ikr.rategkr*ikr.xr*ikr.rkr*(x[0]-var.Ek);

}

// Slowly Activating Potassium Current 
void comp_iks (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.xs1ss = iks.Txs1ss[iV]*d2 + iks.Txs1ss[iV+1]*d1;
	//iks.xs2ss = iks.xs1ss;
	iks.tauxs1 = iks.Ttauxs1[iV]*d2 + iks.Ttauxs1[iV+1]*d1;
	iks.tauxs2 = iks.Ttauxs2[iV]*d2 + iks.Ttauxs2[iV+1]*d1;
	iks.KsCa = 1.0+0.6/(1.0+pow(3.8e-5/x[57],1.4));
	iks.ik = iks.Gks*iks.KsCa*x[47]*x[48]*(x[0]-var.Eks);

}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (double x[])
{
        
	MKL_INT iV=0;   
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ik1.k1ss = ik1.Tk1ss[iV]*d2 + ik1.Tk1ss[iV+1]*d1;
	ik1.tauk1 = ik1.Ttauk1[iV]*d2 + ik1.Ttauk1[iV+1]*d1;
	ik1.rk1 = ik1.Trk1[iV]*d2 + ik1.Trk1[iV+1]*d1;

	ik1.ik = ik1.Gk1*ik1.rategk1*x[49]*ik1.rk1*(x[0]-var.Ek);

}

// Sodium-Calcium Exchanger V-S

void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.hca=var.Thca[iV]*d2 + var.Thca[iV+1]*d1;
	var.hna=var.Thna[iV]*d2 + var.Thna[iV+1]*d1;

	// intracallular space
	ncxi.h1 = 1.0 + (x[53]/var.kna3)*(1.0 + var.hna);
	ncxi.h2 = x[53]*var.hna/(var.kna3*ncxi.h1);
	ncxi.h3 = 1.0/ncxi.h1;
	ncxi.h4 = 1.0 + (x[53]/var.kna1)*(1.0+x[53]/var.kna2);
	ncxi.h5 = x[53]*x[53]/(ncxi.h4*var.kna1*var.kna2);
	ncxi.h6 = 1.0/ncxi.h4;
	ncxi.h7 = 1.0+(var.nao/var.kna3)*(1.0+1.0/var.hna);
	ncxi.h8 = var.nao/(var.kna3*var.hna*ncxi.h7);
	ncxi.h9 = 1.0/ncxi.h7;
	ncxi.h10 = var.kasym+1.0+(var.nao/var.kna1)*(1.0+var.nao/var.kna2);
	ncxi.h11 = var.nao*var.nao/(ncxi.h10*var.kna1*var.kna2);
	ncxi.h12 = 1.0/ncxi.h10;

	ncxi.k1 = ncxi.h12*var.cao*var.kca_on;
	ncxi.k2 = var.kca_off;
	ncxi.k31 = ncxi.h9*var.omega_ca;
	ncxi.k32 = ncxi.h8*var.omega_naca;
	ncxi.k3 = ncxi.k31 + ncxi.k32;
	ncxi.k41 = ncxi.h3*var.omega_ca/var.hca;
	ncxi.k42 = ncxi.h2*var.omega_naca;
	ncxi.k4 = ncxi.k41 + ncxi.k42;
	ncxi.k5 = var.kca_off;
	ncxi.k6 = ncxi.h6*x[57]*var.kca_on;
	ncxi.k7 = ncxi.h5*ncxi.h2*var.omega_na;
	ncxi.k8 = ncxi.h8*ncxi.h11*var.omega_na;
	
	ncxi.x1 = ncxi.k2*ncxi.k4*(ncxi.k7+ncxi.k6)+ncxi.k5*ncxi.k7*(ncxi.k2+ncxi.k3);
	ncxi.x2 = ncxi.k1*ncxi.k7*(ncxi.k4+ncxi.k5)+ncxi.k4*ncxi.k6*(ncxi.k1+ncxi.k8);
	ncxi.x3 = ncxi.k1*ncxi.k3*(ncxi.k7+ncxi.k6)+ncxi.k8*ncxi.k6*(ncxi.k2+ncxi.k3);
	ncxi.x4 = ncxi.k2*ncxi.k8*(ncxi.k4+ncxi.k5)+ncxi.k3*ncxi.k5*(ncxi.k1+ncxi.k8);

	ncxi.E1 = ncxi.x1/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E2 = ncxi.x2/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E3 = ncxi.x3/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E4 = ncxi.x4/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	

	ncxi.allo = 1.0/(1.0+(var.km_ca_act/x[57])*(var.km_ca_act/x[57]));

	ncxi.jnaca_na = 3.0*(ncxi.E4*ncxi.k7 - ncxi.E1*ncxi.k8) + ncxi.E3*ncxi.k42 - ncxi.E2*ncxi.k32;
	ncxi.jnaca_ca = ncxi.E2*ncxi.k2 - ncxi.E1*ncxi.k1;
	var.inaca_i = var.Gnaca*0.8*ncxi.allo*(var.zna*ncxi.jnaca_na + var.zca*ncxi.jnaca_ca);
	
	// subspace
	ncxss.h1 = 1.0 + (x[54]/var.kna3)*(1.0 + var.hna);
	ncxss.h2 = x[54]*var.hna/(var.kna3*ncxss.h1);
	ncxss.h3 = 1.0/ncxss.h1;
	ncxss.h4 = 1.0+x[54]/var.kna1*(1.0+x[54]/var.kna2);
	ncxss.h5 = x[54]*x[54]/(ncxss.h4*var.kna1*var.kna2);
	ncxss.h6 = 1.0/ncxss.h4;
	ncxss.h7 = 1.0+var.nao/var.kna3*(1.0+1.0/var.hna);
	ncxss.h8 = var.nao/(var.kna3*var.hna*ncxss.h7);
	ncxss.h9 = 1.0/ncxss.h7;
	ncxss.h10 = var.kasym+1.0+var.nao/var.kna1*(1.0+var.nao/var.kna2);
	ncxss.h11 = var.nao*var.nao/(ncxss.h10*var.kna1*var.kna2);
	ncxss.h12 = 1.0/ncxss.h10;

	ncxss.k1 = ncxss.h12*var.cao*var.kca_on;
	ncxss.k2 = var.kca_off;
	ncxss.k31 = ncxss.h9*var.omega_ca;
	ncxss.k32 = ncxss.h8*var.omega_naca;
	ncxss.k3 = ncxss.k31 + ncxss.k32;
	ncxss.k41 = ncxss.h3*var.omega_ca/var.hca;
	ncxss.k42 = ncxss.h2*var.omega_naca;
	ncxss.k4 = ncxss.k41 + ncxss.k42;
	ncxss.k5 = var.kca_off;
	ncxss.k6 = ncxss.h6*x[58]*var.kca_on;
	ncxss.k7 = ncxss.h5*ncxss.h2*var.omega_na;
	ncxss.k8 = ncxss.h8*ncxss.h11*var.omega_na;

	ncxss.x1 = ncxss.k2*ncxss.k4*(ncxss.k7+ncxss.k6)+ncxss.k5*ncxss.k7*(ncxss.k2+ncxss.k3);
	ncxss.x2 = ncxss.k1*ncxss.k7*(ncxss.k4+ncxss.k5)+ncxss.k4*ncxss.k6*(ncxss.k1+ncxss.k8);
	ncxss.x3 = ncxss.k1*ncxss.k3*(ncxss.k7+ncxss.k6)+ncxss.k8*ncxss.k6*(ncxss.k2+ncxss.k3);
	ncxss.x4 = ncxss.k2*ncxss.k8*(ncxss.k4+ncxss.k5)+ncxss.k3*ncxss.k5*(ncxss.k1+ncxss.k8);

	ncxss.E1 = ncxss.x1/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E2 = ncxss.x2/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E3 = ncxss.x3/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E4 = ncxss.x4/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	

	ncxss.allo = 1.0/(1.0+(var.km_ca_act/x[58])*(var.km_ca_act/x[58]));

	ncxss.jnaca_na = 3.0*(ncxss.E4*ncxss.k7 - ncxss.E1*ncxss.k8) + ncxss.E3*ncxss.k42 - ncxss.E2*ncxss.k32;
	ncxss.jnaca_ca = ncxss.E2*ncxss.k2 - ncxss.E1*ncxss.k1;
	var.inaca_ss = var.Gnaca*0.2*ncxss.allo*(var.zna*ncxss.jnaca_na + var.zca*ncxss.jnaca_ca);

    var.inaca = var.inaca_i+var.inaca_ss;
    
}

// Sodium-Potassium Pump

void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.knai = inak.Tknai[iV]*d2 + inak.Tknai[iV+1]*d1;
	inak.knao = inak.Tknao[iV]*d2 + inak.Tknao[iV+1]*d1;

	inak.P = inak.SigP/(1.0 + (inak.H/inak.HP) + (x[53]/inak.nap) + (x[55]/inak.kp));

	inak.a1 = inak.kp1*(x[53]/inak.knai)*(x[53]/inak.knai)*(x[53]/inak.knai)/
				((1.0+x[53]/inak.knai)*(1.0+x[53]/inak.knai)*(1.0+x[53]/inak.knai)+(1+x[55]/inak.kki)*(1+x[55]/inak.kki)-1.0);
	inak.b2 = inak.km2*(var.nao/inak.knao)*(var.nao/inak.knao)*(var.nao/inak.knao)/
	 			((1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)+(1.0+var.ko/inak.kko)*(1.0+var.ko/inak.kko)-1.0);
	inak.a3 = inak.kp3*(var.ko/inak.kko)*(var.ko/inak.kko)/
				((1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)+(1.0+var.ko/inak.kko)*(1.0+var.ko/inak.kko)-1.0);
	inak.b3 = inak.km3*inak.P*inak.H/(1.0+inak.MgATP/inak.k_MgATP);
	inak.b4 = inak.km4*(x[55]/inak.kki)*(x[55]/inak.kki)/
				((1.0+x[53]/inak.knai)*(1.0+x[53]/inak.knai)*(1.0+x[53]/inak.knai)+(1.0+x[55]/inak.kki)*(1.0+x[55]/inak.kki)-1.0);

	inak.x1 = inak.a4*inak.a1*inak.a2+inak.b2*inak.b4*inak.b3+inak.a2*inak.b4*inak.b3+inak.b3*inak.a1*inak.a2;
	inak.x2 = inak.b2*inak.b1*inak.b4+inak.a1*inak.a2*inak.a3+inak.a3*inak.b1*inak.b4+inak.a2*inak.a3*inak.b4;
	inak.x3 = inak.a2*inak.a3*inak.a4+inak.b3*inak.b2*inak.b1+inak.b2*inak.b1*inak.a4+inak.a3*inak.a4*inak.b1;
	inak.x4 = inak.b4*inak.b3*inak.b2+inak.a3*inak.a4*inak.a1+inak.b2*inak.a4*inak.a1+inak.b3*inak.b2*inak.a1;

	inak.E1 = inak.x1/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E2 = inak.x2/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E3 = inak.x3/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E4 = inak.x4/(inak.x1+inak.x2+inak.x3+inak.x4);

	inak.jna = 3.0*(inak.E1*inak.a3 - inak.E2*inak.b3);
	inak.jk = 2.0*(inak.E4*inak.b1 - inak.E3*inak.a1);

	inak.inak = inak.G*(var.zna*inak.jna + var.zk*inak.jk);

}

// Sarcolemmal Ca Pump 

void comp_ipca (double x[])
{

	ipca.ca = ipca.G*x[57]/(ipca.km + x[57]);

}

// K Background Current
void comp_ikb (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikb.xkb = ikb.Txkb[iV]*d2 + ikb.Txkb[iV+1]*d1;

	ikb.k = ikb.G*ikb.xkb*(x[0] - var.Ek);
}

// Ca Background Current 

void comp_icab (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	icab.exp = icab.Texp[iV]*d2 + icab.Texp[iV+1]*d1;
	
	icab.ca = icab.pcab*4.0*x[0]*F*F/R/T*(icab.gacai*x[57]*icab.exp-icab.gacao*var.cao)/(icab.exp-1.0);

}

// Na Background Current 

void comp_inab (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inab.exp = inab.Texp[iV]*d2 + inab.Texp[iV+1]*d1;

	inab.na = inab.pnab*x[0]*F*F/R/T*(x[53]*inab.exp-var.nao)/(inab.exp-1.0);

}

void comp_CaMK (double x[])
{
	CaMK.bound = CaMK.z*(1.0-x[50])/(1.0+CaM.Km/x[58]);

	CaMK.active = CaMK.bound + x[50];

}

void comp_diff (double x[])
{

	jdiff.na = (x[54]-x[53])/jdiff.tau_na;
	jdiff.k = (x[56]-x[55])/jdiff.tau_k;
	jdiff.ca = (x[58]-x[57])/jdiff.tau_ca;

}

void comp_jrel (double x[])
{

	jrel.NPss = jrel.p*(-jrel.a*ical.ica/(1.0+pow((1.5/x[60]),8.0)));
	jrel.tau_NP = jrel.b_tau/(1.0+(0.0123/x[60]));
	if (jrel.tau_NP < 0.001){
		jrel.tau_NP = 0.001;
		//printf("jrel tauNP is small\n");
	}

	jrel.CaMKss = jrel.p*(-jrel.a_CaMK*ical.ica/(1.0+pow((1.5/x[60]),8.0)));
	jrel.tau_CaMK = jrel.b_tau_CaMK/(1.0+(0.0123/x[60]));
	if (jrel.tau_CaMK < 0.001){
		jrel.tau_CaMK = 0.001;
	}

	jrel.pCaMK = 1.0/(1.0+CaMK.Km/CaMK.active);

	jrel.ca = (1.0-jrel.pCaMK)*x[51]+jrel.pCaMK*x[52];

}

void comp_jup (double x[])
{

	jup.np = jup.p*(0.004375*x[57]/(0.00092+x[57]));
	jup.CaMK = jup.p*((1.0+jup.dCaMK)*0.004375*x[57]/(0.00092 - jup.dKm_PLB + x[57]));
	
	jup.pCaMK = 1.0/(1.0 + CaMK.Km/CaMK.active);

	jup.leak = 0.0039375*x[59]/15.0;

	jup.ca = (1.0-jup.pCaMK)*jup.np + jup.pCaMK*jup.CaMK - jup.leak;

}

void comp_jtr (double x[])
{

	jtr.ca = (x[59]-x[60])/jtr.tau;

}

void comp_concentration (double x[])
{
	var.b_Ca_i = 1.0/(1.0+var.cmdnbar*var.kmcmdn/((var.kmcmdn+x[57])*(var.kmcmdn+x[57]))+var.trpnbar*var.kmtrpn/((var.kmtrpn+x[57])*(var.kmtrpn+x[57])));
	var.b_Ca_ss = 1.0/(1.0+var.bsrbar*var.kmbsr/((var.kmbsr+x[58])*(var.kmbsr+x[58]))+var.bslbar*var.kmbsl/((var.kmbsl+x[58])*(var.kmbsl+x[58])));
	var.b_Ca_jsr = 1.0/(1.0+var.csqnbar*var.kmcsqn/((var.kmcsqn+x[60])*(var.kmcsqn+x[60])));

}


// Reversal potentials */

void comp_reversal_potential(double x[])
{
	var.Ena = var.RTonF*log(var.nao/x[53]);
	var.Ek = var.RTonF*log(var.ko/x[55]);
	var.Eks = var.RTonF*log((var.ko+var.prnak*var.nao)/(x[55]+var.prnak*x[53]));
	
	//printf("Ena=%lf, Ek=%lf, Eks=%lf\n",var.Ena,var.Ek,var.Eks);
}

