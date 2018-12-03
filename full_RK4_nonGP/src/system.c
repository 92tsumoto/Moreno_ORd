#include "syspara.h"

void function(double x[],double f[],double t)
{

	int i;
	comp_reversal_potential(x);
	comp_CaMK(x);
	comp_ina(x);
	comp_ito(x);
	comp_ical(x);
	comp_ikr(x);
	comp_iks(x);
	comp_ik1(x);
	comp_inaca(x);
	comp_inak(x);
	comp_ipca(x);
	comp_ikb(x);
	comp_icab(x);
	comp_inab(x);
	comp_diff(x);
	comp_jrel(x);
	comp_jup(x);
	comp_jtr(x);
	comp_concentration(x);
	
	var.Ina_i_total = ina.total + inab.na + 3.0*inak.inak + 3.0*var.inaca_i;
	var.Ina_ss_total = ical.icana + 3.0*var.inaca_ss;
	var.Ik_i_total = ito.ik + ikr.ik + iks.ik + ik1.ik + ikb.k - 2.0*inak.inak + var.Istim;
	var.Ik_ss_total = ical.icak;
	var.Ica_i_total = ipca.ca + icab.ca - 2.0*var.inaca_i;
	var.Ica_ss_total = ical.ica - 2.0*var.inaca_ss;
	var.Itotal = ina.total + ito.ik + ical.ica + ical.icana + ical.icak + ikr.ik + iks.ik + ik1.ik 
					+ var.inaca + inak.inak + inab.na + icab.ca + ikb.k + ipca.ca + var.Istim;

	f[0] = -var.Itotal;
	//Fast sodium current from Moreno et al., Sci Transl Med, 2011.
	if(var.drug_type == 0){
	//Drug Free States
		f[1] = ina.a13*x[2] + ina.b2*x[7] + ina.koff*x[9] + ina.k_off*x[18] + ina.bx*x[8] - x[1]* (ina.b13 + ina.a2  + ina.kon + ina.k_on + ina.ax);	// O
		f[2] = ina.a12*x[3] + ina.a3*x[7] + ina.b13*x[1] + ina.kcoff*x[10] + ina.kc_off*x[19] - x[2]*(ina.b12 + ina.b3 + ina.a13  + ina.kcon + ina.kc_on); //C1
		f[3] = ina.a11*x[4] + ina.a3*x[6] + ina.b12*x[2]  + ina.kcoff*x[11] + ina.kc_off*x[20] - x[3]*(ina.b11 + ina.b3 + ina.a12  + ina.kcon + ina.kc_on); //C2
		f[4] = ina.a3*x[5] + ina.b11*x[3]  + ina.kcoff*x[12] + ina.kc_off*x[21] - x[4]*(ina.b3 + ina.a11 + ina.kcon + ina.kc_on); //C3
		f[5] = ina.b3*x[4] + ina.b11*x[6] + ina.ki_off*x[23] - x[5]*(ina.a11 + ina.a3 + ina.ki_on); //IC3
		f[6] = ina.a11*x[5] + ina.b3*x[3] + ina.b12*x[7] + ina.ki_off*x[24] - x[6]*(ina.b11 + ina.a3 + ina.a12 + ina.ki_on); //IC2
		f[7] = ina.a12*x[6] + ina.b3*x[2] + ina.a2*x[1]+ ina.ki_off*x[25] - x[7]*(ina.b12 + ina.b2 + ina.a3 + ina.ki_on); //IF
		f[8] = ina.ax*x[1]+ ina.ki_off*x[22] - x[8]*(ina.bx + ina.ki_on); //IS
	//Charged Drug Bound States
		f[9] = ina.kon*x[1]+ ina.a13c*x[10] + ina.bx1*x[13] + ina.b22*x[16] - x[9]*(ina.koff + ina.b13c + ina.a22 + ina.ax1 ); //DO
		f[10] = ina.kcon*x[2] + ina.a12*x[11] + ina.a33*x[16] + ina.b13c*x[9] - x[10]*(ina.kcoff + ina.b12 + ina.b33 + ina.a13c ); //DC1
		f[11] = ina.kcon*x[3] + ina.a11*x[12] + ina.a33*x[15] + ina.b12*x[10] - x[11]*(ina.kcoff + ina.b11 + ina.b33 + ina.a12 ); //DC2
		f[12] = ina.kcon*x[4] + ina.b11*x[11] + ina.a33*x[14]  - x[12]*(ina.kcoff+ ina.b33 + ina.a11 ); //DC3
		f[13] = ina.ax1*x[9] - ina.bx1*x[13]; //DOS
		f[14] = ina.b33*x[12] + ina.b11*x[15] - x[14]*(ina.a11 + ina.a33); // DIC3
		f[15] = ina.b33*x[11] + ina.a11*x[14] + ina.b12*x[16] - x[15]*(ina.a33 + ina.b11 + ina.a12); //DIC2
		f[16] = ina.b33*x[10] + ina.a12*x[15] + ina.b44*x[17] + ina.a22*x[9] - x[16]*(ina.a33 + ina.b12 + ina.a44 + ina.b22);// DIF
		f[17] = ina.a44*x[16] - x[17]*ina.b44;	//DIM1
	//Neutral Drug Bound States
		f[18] = ina.k_on*x[1]+ ina.a13n*x[19] + ina.b_22*x[25] + ina.bx2*x[22] - x[18]*(ina.k_off + ina.b13n + ina.a_22 + ina.ax2 ); //D_O
		f[19] = ina.kc_on*x[2] + ina.a12*x[20] + ina.a_33*x[25] + ina.b13n*x[18] - x[19]*(ina.kc_off + ina.b12 + ina.b_33 + ina.a13n ); // D_C1
		f[20] = ina.kc_on*x[3] + ina.a11*x[21] + ina.a_33*x[24] + ina.b12*x[19] - x[20]*(ina.kc_off + ina.b11 + ina.b_33 + ina.a12 ); // D_C2
		f[21] = ina.kc_on*x[4] + ina.a_33*x[23] + ina.b11*x[20]  - x[21]*(ina.kc_off + ina.b_33 + ina.a11 ); //D_C3
		f[22] = ina.ax2*x[18] + ina.ki_on*x[8] - x[22]*(ina.bx2 + ina.ki_off); // D_OS
		f[23] = ina.b_33*x[21] + ina.b11*x[24] + ina.ki_on*x[5] - x[23]*(ina.a_33 + ina.a11 + ina.ki_off); //D_IC3
		f[24] = ina.b_33*x[20] + ina.a11*x[23] + ina.b12*x[25] + ina.ki_on*x[6] - x[24]*(ina.a_33 + ina.b11 + ina.a12 + ina.ki_off); //D_IC2
		f[25] = ina.b_33*x[19] + ina.a12*x[24] + ina.b_44*x[26] + ina.a_22*x[18] + ina.ki_on*x[7] - x[25]*(ina.a_33 + ina.a_44 + ina.b_22 + ina.b12 + ina.ki_off); // D_IF
		f[26] = ina.a_44*x[25] - x[26]*ina.b_44; // D_IM1
	} else if(var.drug_type == 1){
	//Drug Free States
		f[1] = ina.a13*x[2] + ina.b2*x[7] + ina.koff*x[9] + ina.k_off*x[18] + ina.bx*x[8] - x[1]* (ina.b13 + ina.a2  + ina.kon + ina.k_on + ina.ax);	
		f[2] = ina.a12*x[3] + ina.a3*x[7] + ina.b13*x[1] + ina.kcoff*x[10] + ina.kc_off*x[19] - x[2]*(ina.b12 + ina.b3 + ina.a13  + ina.kcon + ina.kc_on);
		f[3] = ina.a11*x[4] + ina.a3*x[6] + ina.b12*x[2]  + ina.kcoff*x[11] + ina.kc_off*x[20] - x[3]*(ina.b11 + ina.b3 + ina.a12  + ina.kcon + ina.kc_on);
		f[4] = ina.a3*x[5] + ina.b11*x[3]  + ina.kcoff*x[12] + ina.kc_off*x[21] - x[4]*(ina.b3 + ina.a11 + ina.kcon + ina.kc_on);
		f[5] = ina.b3*x[4] + ina.b11*x[6] + ina.ki_off*x[23] - x[5]*(ina.a11 + ina.a3 + ina.ki_on);
		f[6] = ina.a11*x[5] + ina.b3*x[3] + ina.b12*x[7] + ina.ki_off*x[24] - x[6]*(ina.b11 + ina.a3 + ina.a12 + ina.ki_on);
		f[7] = ina.a12*x[6] + ina.b3*x[2] + ina.a2*x[1]+ ina.ki_off*x[25] - x[7]*(ina.b12 + ina.b2 + ina.a3 + ina.ki_on);
		f[8] = ina.ax*x[1]+ ina.ki_off*x[22] - x[8]*(ina.bx + ina.ki_on);
	//Charged Drug Bound States
		f[9] = ina.kon*x[1]+ ina.a13c*x[10] + ina.bx1*x[13] + ina.b22*x[16] - x[9]*(ina.koff + ina.b13c + ina.a22 + ina.ax1 );
		f[10] = ina.kcon*x[2] + ina.a12*x[11] + ina.a33*x[16] + ina.b13c*x[9] - x[10]*(ina.kcoff + ina.b12 + ina.b33 + ina.a13c );
		f[11] = ina.kcon*x[3] + ina.a11*x[12] + ina.a33*x[15] + ina.b12*x[10] - x[11]*(ina.kcoff + ina.b11 + ina.b33 + ina.a12 );
		f[12] = ina.kcon*x[4] + ina.b11*x[11] + ina.a33*x[14]  - x[12]*(ina.kcoff+ ina.b33 + ina.a11 );
		f[13] = ina.ax1*x[9] - ina.bx1*x[13];
		f[14] = ina.b33*x[12] + ina.b11*x[15] - x[14]*(ina.a11 + ina.a33);
		f[15] = ina.b33*x[11] + ina.a11*x[14] + ina.b12*x[16] - x[15]*(ina.a33 + ina.b11 + ina.a12);
		f[16] = ina.b33*x[10] + ina.a12*x[15] + ina.a22*x[9] - x[16]*(ina.a33 + ina.b12 + ina.b22);
		f[17] = 0.0;
	//Neutral Drug Bound States
		f[18] = ina.k_on*x[1]+ ina.a13n*x[19] + ina.b_22*x[25] + ina.bx2*x[22] - x[18]*(ina.k_off + ina.b13n + ina.a_22 + ina.ax2 );
		f[19] = ina.kc_on*x[2] + ina.a12*x[20] + ina.a_33*x[25] + ina.b13n*x[18] - x[19]*(ina.kc_off + ina.b12 + ina.b_33 + ina.a13n );
		f[20] = ina.kc_on*x[3] + ina.a11*x[21] + ina.a_33*x[24] + ina.b12*x[19] - x[20]*(ina.kc_off + ina.b11 + ina.b_33 + ina.a12 );
		f[21] = ina.kc_on*x[4] + ina.a_33*x[23] + ina.b11*x[20]  - x[21]*(ina.kc_off + ina.b_33 + ina.a11 );
		f[22] = ina.ax2*x[18] + ina.ki_on*x[8] - x[22]*(ina.bx2 + ina.ki_off);
		f[23] = ina.b_33*x[21] + ina.b11*x[24] + ina.ki_on*x[5] - x[23]*(ina.a_33 + ina.a11 + ina.ki_off);
		f[24] = ina.b_33*x[20] + ina.a11*x[23] + ina.b12*x[25] + ina.ki_on*x[6] - x[24]*(ina.a_33 + ina.b11 + ina.a12 + ina.ki_off);
		f[25] = ina.b_33*x[19] + ina.a12*x[24] + ina.a_22*x[18] + ina.ki_on*x[7] - x[25]*(ina.a_33 + ina.b_22 + ina.b12 + ina.ki_off);
		f[26] = 0.0;
	}
	//late sodium current
	f[27] = (ina.mlss - x[27])/ina.tauml; // ml
	f[28] = (ina.hlss - x[28])/ina.tauhl; // hl
	f[29] = (ina.hlCaMKss - x[29])/ina.tauhl_CaMK; // hl
	//Transient outward current
	f[30] = (ito.ass - x[30])/ito.taua;
	f[31] = (ito.iss - x[31])/ito.taui_fast;
	f[32] = (ito.iss - x[32])/ito.taui_slow;
	f[33] = (ito.aCaMKss - x[33])/ito.taua;
	f[34] = (ito.iss - x[34])/ito.taui_CaMK_fast;
	f[35] = (ito.iss - x[35])/ito.taui_CaMK_slow;
	// LTCC
	f[36] = (ical.dss - x[36])/ical.taud;
	f[37] = (ical.fss - x[37])/ical.tauf_fast;
	f[38] = (ical.fss - x[38])/ical.tauf_slow;
	f[39] = (ical.fss - x[39])/ical.taufca_fast;
	f[40] = (ical.fss - x[40])/ical.taufca_slow;
	f[41] = (ical.fss - x[41])/ical.taujca;
	f[42] = (ical.fss - x[42])/ical.tauf_CaMK_fast;
	f[43] = (ical.fss - x[43])/ical.taufca_CaMK_fast;
	f[44] = ical.alpha_n*ical.kp2n - x[44]*ical.km2n;
	// Ikr
	f[45] = (ikr.xrss - x[45])/ikr.tauxr_fast;
	f[46] = (ikr.xrss - x[46])/ikr.tauxr_slow;
	// Iks
	f[47] = (iks.xs1ss - x[47])/iks.tauxs1;
	f[48] = (iks.xs1ss - x[48])/iks.tauxs2;
	// Ik1
	f[49] = (ik1.k1ss - x[49])/ik1.tauk1;
	// CaMK
	f[50] = CaMK.a*CaMK.bound*(CaMK.bound+x[50]) - CaMK.b*x[50];
	// Jrel
	f[51] = (jrel.NPss - x[51])/jrel.tau_NP; 
	f[52] = (jrel.CaMKss - x[52])/jrel.tau_CaMK; 
	// [Na]i
	f[53] = -var.Ina_i_total*var.vr1 + jdiff.na*var.vr2;
	// [Na]ss
	f[54] = -var.Ina_ss_total*var.vr3 - jdiff.na;
	// [K]i
	f[55] = -var.Ik_i_total*var.vr1 + jdiff.k*var.vr2;
	// [K]ss
	f[56] = -var.Ik_ss_total*var.vr3 - jdiff.k;
	// [Ca]i
	f[57] = var.b_Ca_i*(-var.Ica_i_total*var.vr4 - jup.ca*var.vr5 + jdiff.ca*var.vr2);
	// [Ca]ss
	f[58] = var.b_Ca_ss*(-var.Ica_ss_total*var.vr6 +jrel.ca*var.vr7 - jdiff.ca);
	// [Ca]nsr
	f[59] = jup.ca - jtr.ca*var.vr8;
	// [Ca]jsr
	f[60] = var.b_Ca_jsr*(jtr.ca-jrel.ca);

	//printf("NPss=%lf,NPp=%lf\n",jrel.NPss,jrel.CaMKss);
	//for(i=0;i<NN;i++){
	//	printf("x[%d]=%e\n",i,f[i]);
	//}
}
