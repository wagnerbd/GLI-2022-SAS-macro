*************************************************;
*       Macro for Quanjer et al., GLI Spiro     *;
*		Eur Respir J. 2012 Dec 40(6):1324-43 	*;
*		Macro Version 1 - April 7, 2013		    *;
*************************************************;

/* B Wagner modified the 2012 GLI macro code for use				*/
/* with the 2022 race neutral eq from 2022 paper					*/
/* doi:https://doi.org/10.1164/rccm.202205-0963OC 					*/

/* created a lookup table starting with version in paper supplement */
*************************************************;
/*		Instructions for Use					 	

1. Save the macro_code.sas vertical_lookup_table2.sas7bdat to a location
accessable by the sas user.

2. Prepare SAS dataset for calculations:
	Required Fields:
		Name		Level		Description
		id 			numeric 	identifier
		Sex 		numeric 	1 , 2 (Males , Females)
		age 		numeric 	
		height		numeric
		fev1		numeric		if not analyzed value must be .
		fvc			numeric		if not analyzed value must be .
		fev1fvc		numeric		if not analyzed value must be .
		fef2575		numeric		if not analyzed value must be .
		fev075		numeric		if not analyzed value must be .
		fev075fvc	numeric		if not analyzed value must be .


3. 	Skip to Macro Call Code and update with required macro variables		
		Required Fields:
		Name		Description
		data_lib 	The physical location of the saved .sas7bdat to be calculated (no quotations)
		ref_lib		The physical location of the saved vertical_lookup_table2.sas7bdat (no quotations)
		table_name	The name of the saved .sas7bdat table
		out_lib		The physical location to output the  spiro_calculations.sas7bdat table (no quotations)

4.	Run
*/


**********************************************************************************;

* 	Macro Definition; %macro spiro 
				(data_lib= x,
				ref_lib= x,
				table_name= x,
				out_lib= x
				);

*****************************************************************************************;
*	Assign Libraries;

libname datalib "&data_lib";
libname reflib "&ref_lib";
libname outlib "&out_lib";

*****************************************************************************************;

*	Calculations for Fev1;
data fev (keep=id sex age height agebound fev1 FEV1_measurement f);
set datalib.&table_name;
length f $9;


tempage = age- floor(age);
if tempage lt 0.25 then
	agebound = floor(age);
else if tempage lt 0.5 then
	agebound = floor(age)+0.25;
else if tempage lt 0.75 then
	agebound = floor(age) + 0.5;
else
	agebound = floor(age) + 0.75;

FEV1_Measurement=fev1;
f='FEV1';
run;

proc sort data=fev;
by sex f agebound;
run;

proc sort data=reflib.verticle_lookup_table2_2022;
by sex f agebound;
run;

data tempfev;
merge fev reflib.verticle_lookup_table2_2022;
by sex f agebound;
if id ne .;
run;

data fev_calc (keep= id sex age height FEV1_Measurement FEV1_M fev1 FEV1_z_score FEV1_percent_predicted FEV1_LLN l m s);
set work.tempfev;
	l= l0 +(l1-l0)*(age-agebound)/(0.25); /*should be set to 0 for all in GLI 2022*/
	m= m0 +(m1-m0)*(age-agebound)/(0.25);
	s= s0 +(s1-s0)*(age-agebound)/(0.25);
	fev1_M=exp(a0 + (a1*log(height)) + (a2*log(age)) + m);
	fev1_S1=exp(p0+(p1*log(age)) + s);
	fev1_L=q0+(q1*log(age)) + l;
	FEV1_z_score=(((FEV1_Measurement/fev1_M)**fev1_L)-1)/(fev1_L*fev1_S1);
	FEV1_percent_predicted=(FEV1_Measurement/fev1_M)*100;
	FEV1_LLN= fev1_M*(((-1.645*fev1_S1*fev1_L)+1)**(1/fev1_L));
run;


*****************************************************************************************;
*	Calculations for Fvc;

data fvc (keep=id sex age height fvc agebound fvc_measurement f);
set datalib.&table_name;
length f $9;
tempage = age- floor(age);
if tempage lt 0.25 then
	agebound = floor(age);
else if tempage lt 0.5 then
	agebound = floor(age)+0.25;
else if tempage lt 0.75 then
	agebound = floor(age) + 0.5;
else
	agebound = floor(age) + 0.75;

fvc_Measurement=fvc;
f='FVC';
run;

proc sort data=reflib.verticle_lookup_table2_2022;
by sex f agebound;
run;

proc sort data=fvc;
by sex f agebound;
run;

data tempfvc;
merge fvc reflib.verticle_lookup_table2_2022;
by sex f agebound;
if id ne .;
run;

data fvc_calc (keep= id sex age height FVC_Measurement FVC_M FVC_z_score FVC_percent_predicted FVC_LLN);
set work.tempfvc;
	l= l0 +(l1-l0)*(age-agebound)/(0.25);
	m= m0 +(m1-m0)*(age-agebound)/(0.25);
	s= s0 +(s1-s0)*(age-agebound)/(0.25);
	fvc_M=exp(a0+(a1*log(height))+ (a2*log(age)) + m);
	fvc_S1=exp(p0 + (p1*log(age)) + s);
	fvc_L=q0+q1*log(age) + l;
	fvc_z_score=(((fvc_Measurement/fvc_M)**fvc_L)-1)/(fvc_L*fvc_S1);
	fvc_percent_predicted=(fvc_Measurement/fvc_M)*100;
	fvc_LLN= fvc_M*(((-1.645*fvc_S1*fvc_L)+1)**(1/fvc_L));
run;

*****************************************************************************************;
*	Calculations for Fev1fvc;
data fev1fvc (keep=id sex age height fev1fvc agebound fev1fvc_measurement f);
set datalib.&table_name;
length f $9;

tempage = age- floor(age);
if tempage lt 0.25 then
	agebound = floor(age);
else if tempage lt 0.5 then
	agebound = floor(age)+0.25;
else if tempage lt 0.75 then
	agebound = floor(age) + 0.5;
else
	agebound = floor(age) + 0.75;

fev1fvc_Measurement=fev1fvc;
f='FEV1FVC';
run;

proc sort data=fev1fvc;
by sex f agebound;
run;

proc sort data=reflib.verticle_lookup_table2_2022;
by sex f agebound;
run;

data tempfev1fvc;
merge fev1fvc reflib.verticle_lookup_table2_2022;
by sex f agebound;
if id ne .;
run;

data fev1fvc_calc (keep= id sex age height fev1fvc_Measurement fev1fvc_M FEV1FVC_z_score FEV1FVC_LLN FEV1FVC_percent_predicted);
set work.tempfev1fvc;
	l= l0 +(l1-l0)*(age-agebound)/(0.25);
	m= m0 +(m1-m0)*(age-agebound)/(0.25);
	s= s0 +(s1-s0)*(age-agebound)/(0.25);
	fev1fvc_M=exp(a0+(a1*log(height))+(a2*log(age))+m);
	fev1fvc_S1=exp(p0 + (p1*log(age)) + s);
	fev1fvc_L=q0+q1*log(age) + l;
	FEV1FVC_z_score=(((fev1fvc_Measurement/fev1fvc_M)**fev1fvc_L)-1)/(fev1fvc_L*fev1fvc_S1);
	FEV1FVC_LLN= fev1fvc_M*((-1.645*fev1fvc_S1*fev1fvc_L)+1)**(1/fev1fvc_L);
	FEV1FVC_percent_predicted=(fev1fvc_Measurement/fev1fvc_M)*100;

run;

*****************************************************************************************;
*	Calculations for Fef2575;
data fef2575 (keep=id sex age height agebound fef2575_measurement f);
set datalib.&table_name;
length f $9;
tempage = age- floor(age);
if tempage lt 0.25 then
	agebound = floor(age);
else if tempage lt 0.5 then
	agebound = floor(age)+0.25;
else if tempage lt 0.75 then
	agebound = floor(age) + 0.5;
else
	agebound = floor(age) + 0.75;

fef2575_Measurement=fef2575;
f='FEF2575';
run;

proc sort data=fef2575;
by sex f agebound;
run;

proc sort data=reflib.verticle_lookup_table2;
by sex f agebound;
run;

data tempfef2575;
merge fef2575 reflib.verticle_lookup_table2;
by sex f agebound;
if id ne .;
run;

data fef2575_calc (keep= id sex age height FEF2575_Measurement FEF2575_M FEF2575_z_score FEF2575_LLN fef2575_percent_predicted);
set work.tempfef2575;
	l= l0 +(l1-l0)*(age-agebound)/(0.25);
	m= m0 +(m1-m0)*(age-agebound)/(0.25);
	s= s0 +(s1-s0)*(age-agebound)/(0.25);
	fef2575_M=exp(a0+(a1*log(height))+(a2*log(age)) + m);
	fef2575_S1=exp(p0 + (p1*log(age)) + s);
	fef2575_L=q0+q1*log(age) + l;
	fef2575_z_score=(((fef2575_Measurement/fef2575_M)**fef2575_L)-1)/(fef2575_L*fef2575_S1);
	fef2575_LLN= fef2575_M*(((-1.645*fef2575_S1*fef2575_L)+1)**(1/fef2575_L));
	fef2575_percent_predicted=(fef2575_Measurement/fef2575_M)*100;
run;

*****************************************************************************************;
*	Calculations for Fef75;
data fef75 (keep=id sex age height agebound fef75_measurement f);
set datalib.&table_name;
length f $9;

tempage = age- floor(age);
if tempage lt 0.25 then
	agebound = floor(age);
else if tempage lt 0.5 then
	agebound = floor(age)+0.25;
else if tempage lt 0.75 then
	agebound = floor(age) + 0.5;
else
	agebound = floor(age) + 0.75;

fef75_Measurement=fef75;
f='FEF75';
run;

proc sort data=fef75;
by sex f agebound;
run;

proc sort data=reflib.verticle_lookup_table2;
by sex f agebound;
run;

data tempfef75;
merge fef75 reflib.verticle_lookup_table2;
by sex f agebound;
if id ne .;
run;

data fef75_calc (keep= id sex age height FEF75_Measurement FEF75_M FEF75_z_score FEF75_LLN fef75_percent_predicted);
set work.tempfef75;
	l= l0 +(l1-l0)*(age-agebound)/(0.25);
	m= m0 +(m1-m0)*(age-agebound)/(0.25);
	s= s0 +(s1-s0)*(age-agebound)/(0.25);
	fef75_M=exp(a0+(a1*log(height))+(a2*log(age)) + m);
	fef75_S1=exp(p0 + (p1*log(age)) + s);
	fef75_L=q0+q1*log(age) + l;
	fef75_z_score=(((fef75_Measurement/fef75_M)**fef75_L)-1)/(fef75_L*fef75_S1);
	fef75_LLN= fef75_M*(((-1.645*fef75_S1*fef75_L)+1)**(1/fef75_L));
	fef75_percent_predicted=(fef75_Measurement/fef75_M)*100;
run;


/*BW note: none of the following were modified from the original code, ethnicity was not included*/
*****************************************************************************************;
*	Calculations for Fev075;

data fev075 (keep=id sex age height fev075_measurement f);
set datalib.&table_name;
age=round(age, 0.25);
fev075_Measurement=fev075;
f='FEV075';
run;

proc sort data=fev075;
by sex f age;
run;

data fev075_calc (keep= id sex age height FEV075_Measurement FEV075_M FEV075_z_score FEV075_LLN FEV075_percent_predicted);
set work.fev075;
if age >= 3 and age <= 7 then
do;
if sex eq 1 then
	do;
		fev075_M=exp(-9.2947 + 1.9527*log(height)+0.0295*age);
		fev075_S1=exp(- 1.7927 - 0.0666*age);
		fev075_L=1;
		fev075_z_score=(((fev075_Measurement/fev075_M)**fev075_L)-1)/(fev075_L*fev075_S1);
		fev075_LLN= fev075_M*(((-1.645*fev075_S1*fev075_L)+1)**(1/fev075_L));
		fev075_percent_predicted=(fev075_Measurement/fev075_M)*100;
	
	end;
else
	do;
		fev075_M=exp(-10.0239 + 2.1123*log(height)+0.0194*age);
		fev075_S1=exp(-1.6265 + -0.0828*age);
		fev075_L=1;
		fev075_z_score=(((fev075_Measurement/fev075_M)**fev075_L)-1)/(fev075_L*fev075_S1);
		fev075_LLN= fev075_M*(((-1.645*fev075_S1*fev075_L)+1)**(1/fev075_L));
		fev075_percent_predicted=(fev075_Measurement/fev075_M)*100;	
	end;
end;		
run;


*****************************************************************************************;
*	Calculations for Fev075fvc;

data fev075fvc (keep=id sex age height fev075fvc_measurement f);
set datalib.&table_name;
length f $7;
age=round(age, 0.25);
fev075fvc_measurement=fev075fvc;
f='FEV075FVC';
run;

proc sort data=fev075fvc;
by sex f age;
run;

data fev075fvc_calc (keep= id sex age height FEV075fvc_Measurement FEV075fvc_m FEV075fvc_z_score FEV075fvc_LLN FEV075fvc_percent_predicted);
set work.fev075fvc;
if age >= 3 and age <= 7 then
do;
if sex eq 1 then
	do;
		fev075fvc_M=exp(0.5174 + -0.094*log(height)-0.1324*log(age));
		fev075fvc_S1=exp(-2.457);
		fev075fvc_L=2.97;
		fev075fvc_z_score=(((fev075fvc_measurement/fev075fvc_M)**fev075fvc_L)-1)/(fev075fvc_L*fev075fvc_S1);
		fev075fvc_LLN= fev075fvc_M*(((-1.645*fev075fvc_S1*fev075fvc_L)+1)**(1/fev075fvc_L));
		fev075fvc_percent_predicted=(fev075fvc_Measurement/fev075fvc_M)*100;	
	end;
else
	do;
		fev075fvc_M=exp(0.4599 + -0.0891*log(height)+-0.1*log(age));
		fev075fvc_S1=exp(-2.545);
		fev075fvc_L=3.564;
		fev075fvc_z_score=(((fev075fvc_measurement/fev075fvc_M)**fev075fvc_L)-1)/(fev075fvc_L*fev075fvc_S1);
		fev075fvc_LLN= fev075fvc_M*(((-1.645*fev075fvc_S1*fev075fvc_L)+1)**(1/fev075fvc_L));
		fev075fvc_percent_predicted=(fev075fvc_Measurement/fev075fvc_M)*100;	
	end;
end;
run;

*****************************************************************************************;
*	Combine calculated results;

proc sort data=fev_calc;
by id;
run;

proc sort data=fvc_calc;
by id;
run;

proc sort data=fev1fvc_calc;
by id;
run;

proc sort data=fef2575_calc;
by id;
run;

proc sort data=fev075_calc;
by id;
run;

proc sort data=fef75_calc;
by id;
run;

proc sort data=fev075fvc_calc;
by id;
run;

data outlib.spiro_calculations;
merge work.fev_calc work.fvc_calc work.fev1fvc_calc work.fev075_calc work.fev075fvc_calc work.fef2575_calc work.fef75_calc;
by id;
run;

*****************************************************************************************;
*****************************************************************************************;
%mend;

**********************************************************************************;

* 	Macro Call Code;
%spiro	(data_lib= x,
				ref_lib= x,
				table_name= x,
				out_lib= x
		);
