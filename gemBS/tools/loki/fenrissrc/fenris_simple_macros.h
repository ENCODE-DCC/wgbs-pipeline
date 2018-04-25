#define SCALE_RFUNC(XX,PP,ZZ,II,NG) { \
	PP=simple_rf[(XX)]; \
	for(ZZ=0.0,II=0;II<(NG);II++) ZZ+=*PP++; \
	if(ZZ<RESCALE_LIMIT) { \
		ZZ=1.0/ZZ; \
		PP=simple_rf[(XX)]; \
		for(II=0;II<(NG);II++) (*PP++)*=ZZ; \
	} \
}

#define ADD_INFO_KIDS4 { \
	z1=1.0; \
	for(kk=0;kk<n_off;kk++) { \
		p2=kval[kk]; \
		z1*=.25*(p2[g1]+p2[g2]+p2[g3]+p2[g4]); \
	} \
}
	
#define ADD_INFO_KIDS1 { \
	z1=1.0; \
	for(kk=0;kk<n_off;kk++) { \
		p2=kval[kk]; \
		z1*=p2[g1]; \
	} \
}

#define ADD_INFO_KIDS3 { \
	z1=1.0; \
	for(kk=0;kk<n_off;kk++) { \
		p2=kval[kk]; \
		z1*=.25*(p2[g1]+p2[g4])+.5*p2[g2]; \
	} \
}

#define ADD_INFO_KIDS2(xx,yy) { \
	z1=1.0; \
	for(kk=0;kk<n_off;kk++) { \
		p2=kval[kk]; \
		z1*=.5*(p2[xx]+p2[yy]); \
	} \
}

#define HANDLE_SIRE_DAM { \
	z2=z1*z*(*p1++); \
	z3=z1*za*(*p1a++); \
	z4=z2+z3; \
	if(sire_rfp>=0) { \
		p2=simple_rf[sire_rfp]; \
		if(z2) p2[ix1]+=z2/prf[ix1]; \
		if(z3) p2[ix2]+=z3/prf[ix2]; \
	} \
	if(dam_rfp>=0) { \
		p2=simple_rf[dam_rfp]; \
		if(z2) p2[ix2]+=z2/mrf[ix2]; \
		if(z3) p2[ix1]+=z3/mrf[ix1]; \
	} \
	val[ix1]+=z2; \
	val1[ix2]+=z2; \
	val[ix2]+=z3; \
	val1[ix1]+=z3; \
	pp+=z4; \
}

#define HANDLE_SIRE_DAMa { \
	z2=z1*z*(*p1++); \
	if(z2>0.0) { \
		if(sire_rfp>=0) { \
			p2=simple_rf[sire_rfp]; \
			p2[ix1]+=z2/prf[ix1]; \
		} \
		if(dam_rfp>=0) { \
			p2=simple_rf[dam_rfp]; \
			p2[ix1]+=z2/mrf[ix1]; \
		} \
		val[ix1]+=z2; \
		val1[ix1]+=z2; \
		pp+=z2; \
	} \
	z4=z2; \
}

#define HANDLE_KIDS4 { \
	if(z4>0.0) { \
		for(kk=peel_start;kk<n_off1;kk++) { \
			p2=simple_rf[rfp[2+kk]]; \
			p2a=krf[kk]; \
			z5=z4/(p2a[g1]+p2a[g2]+p2a[g3]+p2a[g4]); \
			p2[g1]+=z5; \
			p2[g2]+=z5; \
			p2[g3]+=z5; \
			p2[g4]+=z5; \
		} \
		for(kk=0;kk<n_off;kk++) { \
			p2=kpost[kk]; \
			p2a=kval[kk]; \
			z5=z4/(p2a[g1]+p2a[g2]+p2a[g3]+p2a[g4]); \
			p2[g1]+=z5*p2a[g1]; \
			p2[g2]+=z5*p2a[g2]; \
			p2[g3]+=z5*p2a[g3]; \
			p2[g4]+=z5*p2a[g4]; \
		} \
	} \
} 

#define HANDLE_KIDS1 { \
	if(z4>0.0) { \
		for(kk=peel_start;kk<n_off1;kk++) { \
			p2=simple_rf[rfp[2+kk]]; \
			p2a=krf[kk]; \
			z5=z4/p2a[g1]; \
			p2[g1]+=z5; \
		} \
		for(kk=0;kk<n_off;kk++) { \
			p2=kpost[kk]; \
			p2[g1]+=z4; \
		} \
	} \
} 

#define HANDLE_KIDS3 { \
	if(z4>0.0) { \
		for(kk=peel_start;kk<n_off1;kk++) { \
			p2=simple_rf[rfp[2+kk]]; \
			p2a=krf[kk]; \
			z5=z4/(.25*(p2a[g1]+p2a[g4])+.5*p2a[g2]); \
			p2[g1]+=.25*z5; \
			p2[g2]+=.5*z5; \
			p2[g4]+=.25*z5; \
		} \
		for(kk=0;kk<n_off;kk++) { \
			p2=kpost[kk]; \
			p2a=kval[kk]; \
			z5=z4/(.25*(p2a[g1]+p2a[g4])+.5*p2a[g2]); \
			p2[g1]+=.25*z5*p2a[g1]; \
			p2[g2]+=.5*z5*p2a[g2]; \
			p2[g4]+=.25*z5*p2a[g4]; \
		} \
	} \
} 

#define HANDLE_KIDS2(xx,yy) { \
	if(z4>0.0) { \
		for(kk=peel_start;kk<n_off1;kk++) { \
			p2=simple_rf[rfp[2+kk]]; \
			p2a=krf[kk]; \
			z5=z4/(p2a[xx]+p2a[yy]); \
			p2[xx]+=z5; \
			p2[yy]+=z5; \
		} \
		for(kk=0;kk<n_off;kk++) { \
			p2=kpost[kk]; \
			p2a=kval[kk]; \
			z5=z4/(p2a[xx]+p2a[yy]); \
			p2[xx]+=z5*p2a[xx]; \
			p2[yy]+=z5*p2a[yy]; \
		} \
	} \
} 

