/*----------------------------------------------------------------------------------------------------------
* dofile that calculates predicted 10-year CVD risk based on risk models derived using ERFC 
data, with intention for further recalibration using external estimates.
Jan 2017:	Inital script coded for CHD and Stroke endpoints separately
Feb 2019:	Added supplementary code for composite CVD endpoint
Sep 2019:	Added checking of existence of required variables 
				Renamed linear predictor prefixes to lp_m1_* and lp_m2_*
------------------------------------------------------------------------------------------------------------*/
* check that required variables exist, and if not, generate dummy missing variable
local m1vars = "sex ages smallbin hxdiabbin sbp tchol"
local m2vars = "sex ages smallbin hxdiabbin sbp bmi"
local modelvars: list m1vars | m2vars
foreach var in `modelvars' {
	capture confirm variable `var', exact
	if _rc ~= 0 {
		gen `var' = .
	}
}

gen centchol = tchol - 6
gen censbp = sbp - 120
gen cenages = ages - 60
gen cenbmi = bmi - 25

gen lp_m1_chd=.
gen lp_m1_crbv=.
gen lp_m1_cvdd=.
gen lp_m2_chd=.
gen lp_m2_crbv=.
gen lp_m2_cvdd=.
gen erfc_chdr_m1=.
gen erfc_strr_m1=.
gen erfc_cvdd_m1=.
gen erfc_chdr_m2=.
gen erfc_strr_m2=.
gen erfc_cvdd_m2=.

label variable lp_m1_chd "Linear predictor WHO-CHD Model 1"
label variable lp_m1_crbv "Linear predictor WHO-STR Model 1"
label variable lp_m1_cvdd "Linear predictor WHO-CVD Model 1 (direct)"
label variable lp_m2_chd "Linear predictor WHO-CHD Model 2"
label variable lp_m2_crbv "Linear predictor WHO-STR Model 2"
label variable lp_m2_cvdd "Linear predictor WHO-CVD Model 2 (direct)"
label variable erfc_chdr_m1 "ERFC 10-year risk WHO-CHD Model 1"
label variable erfc_strr_m1 "ERFC 10-year risk WHO-STR Model 1"
label variable erfc_cvdd_m1 "ERFC 10-year risk WHO-CVD Model 1 (direct)"
label variable erfc_chdr_m2 "ERFC 10-year risk WHO-CHD Model 2"
label variable erfc_strr_m2 "ERFC 10-year risk WHO-STR Model 2"
label variable erfc_cvdd_m2 "ERFC 10-year risk WHO-CVD Model 2 (direct)"


*total cholesterol model:
replace lp_m1_chd = 0.0719227*cenages ///
				+0.2284944*centchol ///
				+0.0132183*censbp ///
				+0.6410114*hxdiabbin ///
				+0.5638109*smallbin ///
				-0.0045806*cenages*centchol ///
				-0.0001576*cenages*censbp ///
				-0.0124966*cenages*hxdiabbin ///
				-0.0182545*cenages*smallbin if sex==1
			
replace lp_m1_chd = 0.1020713*cenages ///
				+0.2050377*centchol ///
				+0.015823*censbp ///
				+1.070358*hxdiabbin ///
				+1.053223*smallbin ///
				-0.0051932*cenages*centchol ///
				-0.0001378*cenages*censbp ///
				-0.0234174*cenages*hxdiabbin ///
				-0.0332666*cenages*smallbin if sex==2
	
replace lp_m1_crbv = 0.0986578*cenages ///
				+0.029526*centchol ///
				+0.0222629*censbp ///
				+0.6268712*hxdiabbin ///
				+0.4981217*smallbin ///
				+0.00142*cenages*centchol ///
				-0.0004147*cenages*censbp ///
				-0.026302*cenages*hxdiabbin ///
				-0.0150561*cenages*smallbin if sex==1
						
replace lp_m1_crbv = 0.1056632*cenages ///
				+0.0257782*centchol ///
				+0.0206278*censbp ///
				+0.8581998*hxdiabbin ///
				+0.7443627*smallbin ///
				-0.0021387*cenages*centchol ///
				-0.0004897*cenages*censbp ///
				-0.0209826*cenages*hxdiabbin ///
				-0.0200822*cenages*smallbin if sex==2
				
replace lp_m1_cvdd = 0.0772639*cenages ///
				+0.1799545*centchol ///
				+0.0158976*censbp ///
				+0.6313785*hxdiabbin ///
				+0.5346044*smallbin ///
				-0.0054132*cenages*centchol ///
				-0.000165*cenages*censbp ///
				-0.0173188*cenages*hxdiabbin ///
				-0.0196*cenages*smallbin if sex==1
						
replace lp_m1_cvdd = 0.1020502*cenages ///
				+0.1288215*centchol ///
				+0.0182068*censbp ///
				+0.96491*hxdiabbin ///
				+0.9103094*smallbin ///
				-0.0056806*cenages*centchol ///
				-0.0002959*cenages*censbp ///
				-0.0237451*cenages*hxdiabbin ///
				-0.0289875*cenages*smallbin if sex==2
***********************
			
replace erfc_chdr_m1 = 1-0.9540108^exp(lp_m1_chd) if sex==1
replace erfc_chdr_m1 = 1-0.9894159^exp(lp_m1_chd) if sex==2
replace erfc_strr_m1 = 1-0.9849257^exp(lp_m1_crbv) if sex==1
replace erfc_strr_m1 = 1-0.9890569^exp(lp_m1_crbv) if sex==2
replace erfc_cvdd_m1 = 1-0.9420507^exp(lp_m1_cvdd) if sex==1
replace erfc_cvdd_m1 = 1-0.9791647^exp(lp_m1_cvdd) if sex==2

gen erfc_cvdr_m1=1-(1-erfc_chdr_m1)*(1-erfc_strr_m1)
label variable  erfc_cvdr_m1 "ERFC 10-year risk WHO-CVD Model 1"

*Non-laboratory/BMI model

replace lp_m2_chd = 0.073593*cenages ///
				+0.0337219*cenbmi ///
				+0.0133937*censbp ///
				+0.5954767*smallbin ///
				-0.0010432*cenages*cenbmi ///
				-0.0001837*cenages*censbp ///
				-0.0200831*cenages*smallbin if sex==1
			
replace lp_m2_chd = 0.1049418*cenages ///
				+0.0257616*cenbmi ///
				+0.016726*censbp ///
				+1.093132*smallbin ///
				-0.0006537*cenages*cenbmi ///
				-0.0001966*cenages*censbp ///
				-0.0343739*cenages*smallbin if sex==2
						
replace lp_m2_crbv = 0.097674*cenages ///
				+0.0159518*cenbmi ///
				+0.0227294*censbp ///
				+0.4999862*smallbin ///
				-0.0003516*cenages*cenbmi ///
				-0.0004374*cenages*censbp ///
				-0.0153895*cenages*smallbin  if sex==1
					
replace lp_m2_crbv = 0.1046105*cenages ///
				+0.0036406*cenbmi ///
				+0.0216741*censbp ///
				+0.7399405*smallbin ///
				-0.0000129*cenages*cenbmi ///
				-0.0005311*cenages*censbp ///
				-0.0203997*cenages*smallbin  if sex==2

replace lp_m2_cvdd = 0.079508*cenages ///
				+0.0280687*cenbmi ///
				+0.0162707*censbp ///
				+0.5582418*smallbin ///
				-0.0011959*cenages*cenbmi ///
				-0.0002003*cenages*censbp ///
				-0.0215721*cenages*smallbin  if sex==1
				
replace lp_m2_cvdd = 0.1030354*cenages ///
				+0.0152375*cenbmi ///
				+0.0192623*censbp ///
				+0.9285875*smallbin ///
				-0.0003618*cenages*cenbmi ///
				-0.0003564*cenages*censbp ///
				-0.0296622*cenages*smallbin  if sex==2

replace erfc_chdr_m2 = 1-0.9544258^exp(lp_m2_chd) if sex==1
replace erfc_chdr_m2 = 1-0.9887124^exp(lp_m2_chd) if sex==2
replace erfc_strr_m2 = 1-0.9848260^exp(lp_m2_crbv) if sex==1				
replace erfc_strr_m2 = 1-0.9885706^exp(lp_m2_crbv) if sex==2	
replace erfc_cvdd_m2 = 1-0.9431922^exp(lp_m2_cvdd) if sex==1				
replace erfc_cvdd_m2 = 1-0.9783443^exp(lp_m2_cvdd) if sex==2	

gen erfc_cvdr_m2=1-(1-erfc_chdr_m2)*(1-erfc_strr_m2)
label variable  erfc_cvdr_m2 "ERFC 10-year risk WHO-CVD Model 2"

* SK clone/label easy to use risk variables
/* 
clonevar erfc_chdr_m1 = erfc_chdr
clonevar erfc_chdr_m2 = erfcbmi_chdr
clonevar erfc_strr_m1 = erfc_crbvr
clonevar erfc_strr_m2 = erfcbmi_crbvr
clonevar erfc_cvdr_m1 = erfc_cvdr
clonevar erfc_cvdr_m2 = erfcbmi_cvdr
*/
forvalues j = 1/2 {
	label variable erfc_chdr_m`j' "WHO-CHD Model `j'"
	label variable erfc_strr_m`j' "WHO-STR Model `j'"
	label variable erfc_cvdr_m`j' "WHO-CVD Model `j'"
	label variable erfc_cvdd_m`j' "WHO-CVD Model `j'"
}
