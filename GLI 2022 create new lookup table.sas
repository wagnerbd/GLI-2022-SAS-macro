/*create a new lookup table for GLI 2022 race neutral equations */
/* start with supplemental table in paper doi:https://doi.org/10.1164/rccm.202205-0963OC 					*/

/*read in supplemental table*/
/*reformat to match 2012 look-up table*/
/*save as: reflib.verticle_lookup_table2_2022*/

libname reflib "C:\Users\wagnerbd\Documents\growth curve macros";


%macro lookup (sheet, sex, pft, a0, a1, a2, p0, p1, q0, q1, out);
PROC IMPORT OUT= pft 
            DATAFILE= "C:\Users\wagnerbd\Documents\growth curve macros\gli_global_lookuptables_dec6.xlsx" 
            DBMS=xlsx REPLACE;
     GETNAMES=YES;
	 Sheet=&sheet.;
RUN;
proc sort data=pft;
	by descending age;
run;
data &out.;
	set pft;
		length f $9;
		format a0 a1 a2 p0 p1 q0 q1 l0 l1 m0 m1 s0 s1 10.8;
		sex=&sex.;
		f=&pft.;
		agebound=age;
		a0=&a0.;		/*intercept term for M equation*/
		a1=&a1.;		/*parameter for ln(height) in M equation*/
		a2=&a2.;		/*parameter for ln(age) in M equation*/

		p0=&p0.;		/*intercept term for S equation*/
		p1=&p1.;		/*parameter for ln(age) in S equation*/

		q0=&q0.;		/*intercept term for L equation*/
		q1=&q1.;		/*parameter for ln(age) in L equation*/

		m0=M_Spline;
		m1=lag(M_Spline);

		s0=S_Spline;
		s1=lag(S_Spline);

		l0=Lspline;
		l1=lag(Lspline);

		if agebound=95 then do;
			m1=m0; s1=s0; l1=l0;
		end;

	drop E G age M_Spline S_Spline Lspline ;
run;
proc sort data=&out.;
	by age;
run;
data pft; run;
%mend;

%lookup ("Male FEV1", 1, 'FEV1', -11.399108, 2.462664, -0.011394, -2.256278, 0.080729, 1.22703, 0, male_fev);
%lookup ("Male FVC", 1, 'FVC', -12.629131, 2.727421, 0.009174, -2.195595, 0.068466, 0.9346, 0, male_fvc);
%lookup ("Male FEV1 FVC", 1, 'FEV1FVC', 1.022608, -0.218592, -0.027586, -2.882025, 0.068889, 3.8243, -0.3328, male_fev_fvc);

%lookup ("Female FEV1", 2, 'FEV1', -10.901689, 2.385928, -0.076386, -2.364047, 0.129402, 1.21388, 0, female_fev);
%lookup ("Female FVC", 2, 'FVC', -12.055901, 2.621579, -0.035975, -2.310148, 0.120428, 0.899, 0, female_fvc);
%lookup ("Female FEV1 FVC", 2, 'FEV1FVC', 0.9189568, -0.1840671, -0.0461306, -3.171582, 0.144358, 6.6490, -0.9920, female_fev_fvc);


data reflib.verticle_lookup_table2_2022;
	retain sex f a0 a1 a2 p0 p1 q0 q1 agebound l0 l1 m0 m1 s0 s1;
	set male_fev male_fvc male_fev_fvc female_fev female_fvc female_fev_fvc;
run;
