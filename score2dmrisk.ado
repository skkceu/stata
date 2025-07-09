*! Version : 0.00 12Apr2022
*! Authors : Dr Stephen Kaptoge
*!				 Dr Lisa Pennells
*!				 Dr Steph Read
*! Address : Cardiovascular Epidemiology Unit
*!           Dept of Public Health and Primary Care
*!           University of Cambridge
*!           Strangeways Research Laboratory
*!           Worts Causeway
*!           Cambridge CB1 8RN
*!           UK

/*----------------------------------------------------------------------------------------------------------
S Kaptoge: Apr 2022
* Stata program for calculation of SCORE2-DM CVD risk score, with country/region specific recalibration.

ERFC/UKB/CPRD/SCIDIAB-derived SCORE2-DM risk prediction model for estimation of 10-year CVD risk,
with additional diabetes variables and adjustment for non-CVD death as competing risk using IPW weighted 
stratified Cox regression, and reclaibration to four ESC risk regions 

* required variables: 
tchol = total cholesterol in mmol/l
sbp = systolic blood pressure in mmHg
ages = baseline age in years
sex = sex (1 = male, 2 = female)
smallbin = smoking as a binary variable (current = 1, other = 0)
hxdiabbin = history of diabetes as a binary variable (definite = 1, other = 0)
agediab = age at diagnosis of diabetes in years
hba1c_ifcc = glycated haemoglobin (HbA1c) in mmol/mol
egfr_ckdepi = estimated glomerular filtration rate (eGFR) using CKDEpi formula in ml/min/1.73m^2

riskreg =  the risk region (1 = low, 2 = Moderate, 3 = high, 4 = very high) 
that applies to the country of residence of the individuals in the dataset 
(see separate table for risk region classification of each country)

relevant descriptive notes.
Jan 2022:	Inital script coded as stata dofile for risk score calculation in validation data
Apr 2022:	Above converted to current wrapper program for general use
Apr 2022:	Added optional calculation of risk age
Apr 2022:	Added option multdiab(integer 1) to optionally include/exclude consideration of diabetes status
Apr 2022:	Added option nobinary to treat smoking and diabetes as continuous variables if using prevalence
Apr 2022:	Added [if] option to be applied at the end of calculations
------------------------------------------------------------------------------------------------------------*/

capture program drop score2dmrisk
program define score2dmrisk, rclass
	
	version 13
	syntax [if] [, CCODEs riskreg(numlist integer max=1 >=1 <=104) replace percent shortlab rounded(integer 0) ///
	dppscale(real 0) noWEB WEBVERSion(integer 2) DETAILs RISKAGE MULTDIAB(integer 1) MULTEXTRA(integer 1) ///
	FULLADJ(integer 1) nobinary]
	
	marksample touse, novarlist
	
	* preserve data in case calculation fails
	preserve
	* web settings
	local webopt = cond("`web'" == "noweb", 0, `webversion')
	local quinoi = "qui"
	`quinoi' {		/* start of quinoi */
		* drop variables from previous calculations if they exist
		capture ds _SCORE2DM_START - _SCORE2DM_END
		if _rc == 0 {
			drop `r(varlist)'
		}
		
		* check if required predictor variables exist
		local expvars = "sex ages sbp tchol hdl smallbin hxdiabbin agediab hba1c_ifcc lnegfr_ckdepi"
		noi confirm numeric variable `expvars'
		capture confirm variable riskreg, exact
		if _rc~=0 {
			if missing("`riskreg'") {
				noi di _newline as text "Please enter SCORE2-DM risk region (1 = Low, 2 = Moderate, 3 = High, 4 = Very high)"
				qui di _request(_answer)
				local riskreg = `answer'
			}
			local chkriskreg = inlist(`riskreg', 1, 2, 3, 4) | inlist(`riskreg', 101, 102, 103, 104)
			if `chkriskreg' == 0 {
				numlist "1(1)4 101(1)104"
				di as error "Invalid riskreg(`riskreg') specification, enter value among: `r(numlist)'"
				exit 197
			}
			gen riskreg = `riskreg'
			label define riskreg 1 "Europe Low" 2 "Europe Mod" 3 "Europe High" 4 "Europe vHigh", modify
			label define riskreg 101 "Asia Low" 102 "Asia Mod" 103 "Asia High" 104 "Asia vHigh", modify	
			label define riskreg 1 "Europe Low risk region" 2 "Europe Moderate risk region" 3 "Europe High risk region" 4 "Europe Very high risk region", modify
			label define riskreg 101 "Asia Low risk region" 102 "Asia Moderate risk region" 103 "Asia High risk region" 104 "Asia Very high risk region", modify
			label values riskreg riskreg
			label variable riskreg "SCORE2 risk region"
			local orderreg = "order riskreg, after(_SCORE2DM_START)"
		}
		else {
			capture confirm numeric variable riskreg, exact
			if _rc == 0 {
				 if "`replace'" == "replace" & !missing("`riskreg'") { 
					replace riskreg = `riskreg'
				}
			}
			else {
				di as error "variable riskreg is not numeric"
				des riskreg
				exit 197
			}
		}
				
		* generate marker for start of SCORE2-DM added variables
		gen _SCORE2DM_START = .
		`orderreg'		/* order risk region variable if generated above */
		
		* drop existing derived variables
		local dropvars1 = "cages csbp cbmi ctchol chdl cagediab chba1c_ifcc clnegfr_ckdepi clnegfr"
		local dropvars2 = "lp_m1_ep_cv s010_m1_ep_cv scale1 scale2 cal0_m1_risk cal1_m1_risk"
		local dropvars3 = "score2dm"
		local dropvars = "`dropvars1' `dropvars2' `dropvars3'"
		local varabbrev = c(varabbrev)
		nobreak {
			set varabbrev off
			foreach var in `dropvars' {
				capture drop `var'
			}
			set varabbrev `varabbrev'
		}
		* center variables and scale as appropriate
		local eqn_cages = "(ages - 60)/5"
		local eqn_csbp = "(sbp - 120)/20"
		local eqn_ctchol = "(tchol - 6)/1"
		local eqn_chdl = "(hdl - 1.3)/0.5"
		local eqn_cagediab = "(agediab - 50)/5"
		local eqn_chba1c_ifcc = "(hba1c_ifcc - 31)/9.34"
		local eqn_clnegfr_ckdepi = "(lnegfr_ckdepi - 4.5)/0.15"
		local expvars = "ages sbp tchol hdl agediab hba1c_ifcc lnegfr_ckdepi"
		foreach expvar of varlist `expvars' {
			local expvarlab: variable label `expvar'
			local cvar = "c`expvar'"
			local cvareqn = "`eqn_`cvar''"
			gen `cvar' = `cvareqn'
			label variable `cvar' "Centered `expvarlab': `cvareqn'"
		}
		clonevar clnegfr = clnegfr_ckdepi
		local hxdiabbin = cond(inlist("`binary'", "nobinary"), "c.hxdiabbin", "1.hxdiabbin")
		local smallbin = cond(inlist("`binary'", "nobinary"), "c.smallbin", "1.smallbin")
		
		* generate linear predictor and baseline survival from CVD Cox FGR, 
		*sex = 1 (Men)
		generat lp_m1_ep_cv =  + 0.3742*(cages) + 0.6012*(`smallbin') + 0.2777*(csbp) + 0.6457*(`hxdiabbin')*`multdiab' + 0.1458*(ctchol) + -0.2698*(chdl) + -0.0755*(`smallbin'#c.cages) + -0.0255*(c.cages#c.csbp) + -0.0983*(`hxdiabbin'#c.cages)*`multdiab' + -0.0281*(c.cages#c.ctchol) + 0.0426*(c.cages#c.chdl) if sex == 1
		replace lp_m1_ep_cv = lp_m1_ep_cv + -0.0998*(`hxdiabbin'#c.cagediab)*`multdiab' + 0.0955*(chba1c_ifcc)*`multextra' + -0.0591*(clnegfr)*`multextra' + 0.0058*(c.clnegfr#c.clnegfr)*`multextra' + -0.0134*(c.cages#c.chba1c_ifcc)*`multextra' + 0.0115*(c.cages#c.clnegfr)*`multextra' if sex == 1
		if `fulladj' == 1 {
			replace lp_m1_ep_cv = lp_m1_ep_cv + 0.1626*(cages) + -0.1454*(csbp) + -0.0356*(ctchol) + 0.1611*(chdl) + -0.1238*(`smallbin') + -0.0013*(c.cages#c.csbp) + 0.0100*(c.cages#c.ctchol) + -0.0331*(c.cages#c.chdl) + 0.0083*(`smallbin'#c.cages) if sex == 1
		}
		else {
			replace lp_m1_ep_cv = lp_m1_ep_cv + 0.1626*(cages) if sex == 1
		}
		generat s010_m1_ep_cv = 0.9605 if sex == 1
		
		*sex = 2 (Women)
		replace lp_m1_ep_cv =  + 0.4648*(cages) + 0.7744*(`smallbin') + 0.3131*(csbp) + 0.8096*(`hxdiabbin')*`multdiab' + 0.1002*(ctchol) + -0.2606*(chdl) + -0.1088*(`smallbin'#c.cages) + -0.0277*(c.cages#c.csbp) + -0.1272*(`hxdiabbin'#c.cages)*`multdiab' + -0.0226*(c.cages#c.ctchol) + 0.0613*(c.cages#c.chdl) if sex == 2
		replace lp_m1_ep_cv = lp_m1_ep_cv + -0.1180*(`hxdiabbin'#c.cagediab)*`multdiab' + 0.1173*(chba1c_ifcc)*`multextra' + -0.0640*(clnegfr)*`multextra' + 0.0062*(c.clnegfr#c.clnegfr)*`multextra' + -0.0196*(c.cages#c.chba1c_ifcc)*`multextra' + 0.0169*(c.cages#c.clnegfr)*`multextra' if sex == 2
		if `fulladj' == 1 {
			replace lp_m1_ep_cv = lp_m1_ep_cv + 0.1976*(cages) + -0.1709*(csbp) + 0.0126*(ctchol) + 0.1038*(chdl) + -0.1605*(`smallbin') + 0.0111*(c.cages#c.csbp) + 0.0026*(c.cages#c.ctchol) + -0.0426*(c.cages#c.chdl) + -0.0033*(`smallbin'#c.cages) if sex == 2
		}
		else {
			replace lp_m1_ep_cv = lp_m1_ep_cv + 0.1976*(cages) if sex == 2
		}
		replace s010_m1_ep_cv = 0.9776 if sex == 2
		
		* replace with rounded coefficients as published for sensitivity comparisons 
		if `rounded' == 1 {
			*sex = 1 (Men)
			if `fulladj' == 1 {
				replace lp_m1_ep_cv =  ln(1.71)*(cages) + ln(1.61)*(`smallbin') + ln(1.14)*(csbp) + ln(1.91)*(`hxdiabbin')*`multdiab' + ln(1.12)*(ctchol) + ln(0.90)*(chdl) + ln(0.94)*(`smallbin'#c.cages) + ln(0.97)*(c.cages#c.csbp) + ln(0.91)*(`hxdiabbin'#c.cages)*`multdiab' + ln(0.98)*(c.cages#c.ctchol) + ln(1.01)*(c.cages#c.chdl) if sex == 1
				replace lp_m1_ep_cv = lp_m1_ep_cv + ln(0.90)*(`hxdiabbin'#c.cagediab)*`multdiab' + ln(1.10)*(chba1c_ifcc)*`multextra' + ln(0.94)*(clnegfr)*`multextra' + ln(1.01)*(c.clnegfr#c.clnegfr)*`multextra' + ln(0.99)*(c.cages#c.chba1c_ifcc)*`multextra' + ln(1.01)*(c.cages#c.clnegfr)*`multextra' if sex == 1
			}
			else {
				replace lp_m1_ep_cv =  ln(1.45)*(cages) + ln(1.82)*(`smallbin') + ln(1.32)*(csbp) + ln(1.91)*(`hxdiabbin')*`multdiab' + ln(1.16)*(ctchol) + ln(0.76)*(chdl) + ln(0.93)*(`smallbin'#c.cages) + ln(0.98)*(c.cages#c.csbp) + ln(0.91)*(`hxdiabbin'#c.cages)*`multdiab' + ln(0.97)*(c.cages#c.ctchol) + ln(1.04)*(c.cages#c.chdl) if sex == 1
				replace lp_m1_ep_cv = lp_m1_ep_cv + ln(0.90)*(`hxdiabbin'#c.cagediab)*`multdiab' + ln(1.10)*(chba1c_ifcc)*`multextra' + ln(0.94)*(clnegfr)*`multextra' + ln(1.01)*(c.clnegfr#c.clnegfr)*`multextra' + ln(0.99)*(c.cages#c.chba1c_ifcc)*`multextra' + ln(1.01)*(c.cages#c.clnegfr)*`multextra' if sex == 1
				replace lp_m1_ep_cv = lp_m1_ep_cv + ln(1.18)*(cages) if sex == 1
			}
			replace s010_m1_ep_cv = 0.9605 if sex == 1
			
			*sex = 2 (Women)
			if `fulladj' == 1 {
				replace lp_m1_ep_cv =  ln(1.94)*(cages) + ln(1.85)*(`smallbin') + ln(1.15)*(csbp) + ln(2.25)*(`hxdiabbin')*`multdiab' + ln(1.12)*(ctchol) + ln(0.85)*(chdl) + ln(0.89)*(`smallbin'#c.cages) + ln(0.98)*(c.cages#c.csbp) + ln(0.88)*(`hxdiabbin'#c.cages)*`multdiab' + ln(0.98)*(c.cages#c.ctchol) + ln(1.02)*(c.cages#c.chdl) if sex == 2
				replace lp_m1_ep_cv = lp_m1_ep_cv + ln(0.89)*(`hxdiabbin'#c.cagediab)*`multdiab' + ln(1.12)*(chba1c_ifcc)*`multextra' + ln(0.94)*(clnegfr)*`multextra' + ln(1.01)*(c.clnegfr#c.clnegfr)*`multextra' + ln(0.98)*(c.cages#c.chba1c_ifcc)*`multextra' + ln(1.02)*(c.cages#c.clnegfr)*`multextra' if sex == 2
			}
			else {
				replace lp_m1_ep_cv =  ln(1.59)*(cages) + ln(2.17)*(`smallbin') + ln(1.37)*(csbp) + ln(2.25)*(`hxdiabbin')*`multdiab' + ln(1.11)*(ctchol) + ln(0.77)*(chdl) + ln(0.90)*(`smallbin'#c.cages) + ln(0.97)*(c.cages#c.csbp) + ln(0.88)*(`hxdiabbin'#c.cages)*`multdiab' + ln(0.98)*(c.cages#c.ctchol) + ln(1.06)*(c.cages#c.chdl) if sex == 2
				replace lp_m1_ep_cv = lp_m1_ep_cv + ln(0.89)*(`hxdiabbin'#c.cagediab)*`multdiab' + ln(1.12)*(chba1c_ifcc)*`multextra' + ln(0.94)*(clnegfr)*`multextra' + ln(1.01)*(c.clnegfr#c.clnegfr)*`multextra' + ln(0.98)*(c.cages#c.chba1c_ifcc)*`multextra' + ln(1.02)*(c.cages#c.clnegfr)*`multextra' if sex == 2
				replace lp_m1_ep_cv = lp_m1_ep_cv + ln(1.22)*(cages) if sex == 2
			}
			replace s010_m1_ep_cv = 0.9776 if sex == 2
		}
		
		* Calculate uncalibrated 10-year risk calculation from linear predictors and baseline survival
		gen cal0_m1_risk = 1 - s010_m1_ep_cv^exp(lp_m1_ep_cv)
		
		* generate recalibration factors depending on sex and risk region of the individuals
		gen scale1 = .
		gen scale2 = .

		**** Europe: Published work ****
		* Europe low risk region
		replace scale1 = -0.569857 if riskreg==1 & sex==1
		replace scale2 =  0.747593 if riskreg==1 & sex==1
		replace scale1 = -0.737966 if riskreg==1 & sex==2
		replace scale2 =  0.701868 if riskreg==1 & sex==2
		
		* Europe moderate risk region
		replace scale1 = -0.156481 if riskreg==2 & sex==1
		replace scale2 = 0.800869 if riskreg==2 & sex==1
		replace scale1 = -0.314269 if riskreg==2 & sex==2
		replace scale2 = 0.770149 if riskreg==2 & sex==2
		
		* Europe high risk region
		replace scale1 = 0.320717 if riskreg==3 & sex==1
		replace scale2 = 0.936004 if riskreg==3 & sex==1
		replace scale1 = 0.570963 if riskreg==3 & sex==2
		replace scale2 = 0.936893 if riskreg==3 & sex==2
		
		* Europe very high risk region
		replace scale1 = 0.583578 if riskreg==4 & sex==1
		replace scale2 = 0.829368 if riskreg==4 & sex==1
		replace scale1 = 0.941177 if riskreg==4 & sex==2
		replace scale2 = 0.832923 if riskreg==4 & sex==2
		
		* replace with rounded coefficients as published for sensitivity comparisons 
		if `rounded' == 1 | `dppscale' > 0 {
			local dpprec = cond(`rounded' == 1 & `dppscale' == 0, 0.01, `dppscale')
			replace scale1 = round(scale1, `dpprec')
			replace scale2 = round(scale2, `dpprec')
		}
		
		* recalibrate the risk estimates for all individuals:
		gen cal1_m1_risk = 1-exp(-exp(scale1 + scale2*ln(-ln(1 - cal0_m1_risk))))
		clonevar score2dm = cal1_m1_risk
		
		* label variables informatively
		label variable lp_m1_ep_cv "Linear predictor CVD"
		label variable s010_m1_ep_cv "COHORT baseline survior at 10 years"
		label variable cal0_m1_risk "COHORT 10-yr CVD risk"
		label variable scale1 "SCORE2-DM rescaling factor 1 (intercept)"
		label variable scale2 "SCORE2-DM rescaling factor 2 (slope)"
		label variable cal1_m1_risk "SCORE2-DM 10-yr CVD risk"
		label variable score2dm "SCORE2-DM 10-yr CVD risk"
		if "`shortlab'" == "shortlab" {
			label variable lp_m1_ep_cv "LP-CVD"
			label variable s010_m1_ep_cv "UTR-s010"
			label variable cal0_m1_risk "UTR-CVD"
			label variable scale1 "SCORE2-DM-b0"
			label variable scale2 "SCORE2-DM-b1"
			label variable cal1_m1_risk "SCORE2-DM-CVD"
			label variable score2dm "SCORE2-DM-CVD"
		}
		
		if "`percent'" == "percent" {
			local riskvars = "cal?_m?_risk score2dm"
			foreach var of varlist `riskvars' {
				local varlab: variable label `var'
				replace `var' = 100*`var'		/* convert to percent risk */
				label variable `var' "`varlab' (%)"
			}
		}
		
		* calculate risk age if requested
		if !missing("`riskage'") {
			score2dmriskage if !missing(score2dm)		/* call wrapper program below to calculate risk age */
		}

		* generate marker for end of SCORE2 added variables
		gen _SCORE2DM_END = .
	}		/* end of quinoi */
	
	* cancel restore if successful
	restore, not
	
	* keep minimal variables in dataset if nodetails
	if inlist(`webopt', 1, 2) {
		qui ds _SCORE2DM_START - _SCORE2DM_END
		local dropvars = r(varlist)
		local calvars = cond(inlist(`webopt', 1), "cal?_m?_risk score2dm", "score2dm")
		if !missing("`riskage'") {
			local calvars = "`calvars' ages riskage?"
		}
		if "`details'" ~= "details" {
			local keepvars = "_SCORE2DM_START riskreg `calvars' _SCORE2DM_END"
			qui ds `keepvars'
			local keepvars = r(varlist)
			local dropvars: list dropvars - keepvars
			drop `dropvars'
		}
	}
	
	* keep minimum riskage variables if nodetails
	if !missing("`riskage'") {
		if "`details'" ~= "details" {
			foreach k of numlist 2(1)6 {
				capture drop *riskage`k'
			}
		}
	}
	
	* describe and summarize risk variables
	qui recode _SCORE2DM_START - _SCORE2DM_END (* = .) if !(`touse')
	des _SCORE2DM_START - _SCORE2DM_END
	preserve
	qui keep if !missing(score2)
	local calvars = cond(inlist(`webopt', 0, 1), "cal?_m?_risk score2dm", "score2dm")
	if !missing("`riskage'") {
		local calvars = "`calvars' ages riskage?"
	}
	bysort riskreg: tabstat `calvars', by(sex) stats(n mean sd min p25 p50 p75 max) ///
	col(stats) longstub varw(16)
	if !missing("`riskage'") {
		notes `calvars'
	}
	restore
	
end

* wrapper program to calculate risk age for SCORE2-DM
capture program drop score2dmriskage
program define score2dmriskage, rclass
	
	version 13
	syntax [if] [in] [, *]
	marksample touse, novarlist
	
	local dropvars = "riskage1 riskage2 riskage3 riskage4 riskage5 riskage6"
	foreach var in `dropvars' {
		capture drop `var'
	}
	foreach i of numlist 1(1)6 {
		gen riskage`i' = .
		label variable riskage`i' "SCORE2-DM risk age`i'"
	}
	notes riskage1: "SCORE2-DM risk age1: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 120, nonHDL = 3.9 (TCHOL = 5.0 & HDL = 1.1))"
	notes riskage2: "SCORE2-DM risk age2: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.9 (TCHOL = 5.0 & HDL = 1.1))"
	notes riskage3: "SCORE2-DM risk age3: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 140, nonHDL = 3.9 (TCHOL = 5.0 & HDL = 1.1))"
	notes riskage4: "SCORE2-DM risk age4: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 110, nonHDL = 3.4 (TCHOL = 4.5 & HDL = 1.1))"
	notes riskage5: "SCORE2-DM risk age5: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 120, nonHDL = 3.9 (TCHOL = 5.0 & HDL = 1.1))"
	notes riskage6: "SCORE2-DM risk age6: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 4.4 (TCHOL = 5.5 & HDL = 1.1))"
	tempvar score2dmpct
	summ score2dm
	gen `score2dmpct' = cond(r(max) < 1, score2dm*100, score2dm)
	* risk age estimation equations derived separately (see score2dmrisk_test.do)
	replace riskage1 = min(., 44.3150 + 12.3421*ln(`score2dmpct') + 0.0838*`score2dmpct'^1) if riskreg == 1 & sex == 1
	replace riskage1 = min(., 42.4866 + 11.7646*ln(`score2dmpct') + 0.0887*`score2dmpct'^1) if riskreg == 2 & sex == 1
	replace riskage1 = min(., 37.2353 + 12.9862*ln(`score2dmpct') + 0.0013*`score2dmpct'^2) if riskreg == 3 & sex == 1
	replace riskage1 = min(., -0.7265 + 21.5636*ln(`score2dmpct') + 0.0019*`score2dmpct'^2) if riskreg == 4 & sex == 1
	replace riskage1 = min(., 50.2173 + 11.2004*ln(`score2dmpct') + 0.0741*`score2dmpct'^1) if riskreg == 1 & sex == 2
	replace riskage1 = min(., 49.7558 + 10.2392*ln(`score2dmpct') + 0.0764*`score2dmpct'^1) if riskreg == 2 & sex == 2
	replace riskage1 = min(., 44.6856 + 10.4804*ln(`score2dmpct') + 0.0011*`score2dmpct'^2) if riskreg == 3 & sex == 2
	replace riskage1 = min(., 14.7232 + 16.6904*ln(`score2dmpct') + 0.0016*`score2dmpct'^2) if riskreg == 4 & sex == 2
	replace riskage2 = min(., 40.6356 + 13.2785*ln(`score2dmpct') + 0.0899*`score2dmpct'^1) if riskreg == 1 & sex == 1
	replace riskage2 = min(., 38.6747 + 12.6537*ln(`score2dmpct') + 0.0951*`score2dmpct'^1) if riskreg == 2 & sex == 1
	replace riskage2 = min(., 32.9752 + 13.9975*ln(`score2dmpct') + 0.0014*`score2dmpct'^2) if riskreg == 3 & sex == 1
	replace riskage2 = min(., -7.9871 + 23.2673*ln(`score2dmpct') + 0.0020*`score2dmpct'^2) if riskreg == 4 & sex == 1
	replace riskage2 = min(., 47.8112 + 11.7850*ln(`score2dmpct') + 0.0780*`score2dmpct'^1) if riskreg == 1 & sex == 2
	replace riskage2 = min(., 47.3275 + 10.7716*ln(`score2dmpct') + 0.0804*`score2dmpct'^1) if riskreg == 2 & sex == 2
	replace riskage2 = min(., 41.9772 + 11.0400*ln(`score2dmpct') + 0.0011*`score2dmpct'^2) if riskreg == 3 & sex == 2
	replace riskage2 = min(., 10.3927 + 17.5947*ln(`score2dmpct') + 0.0016*`score2dmpct'^2) if riskreg == 4 & sex == 2
	replace riskage3 = min(., 36.3505 + 14.3681*ln(`score2dmpct') + 0.0971*`score2dmpct'^1) if riskreg == 1 & sex == 1
	replace riskage3 = min(., 34.2369 + 13.6873*ln(`score2dmpct') + 0.1026*`score2dmpct'^1) if riskreg == 2 & sex == 1
	replace riskage3 = min(., 28.0013 + 15.1800*ln(`score2dmpct') + 0.0015*`score2dmpct'^2) if riskreg == 3 & sex == 1
	replace riskage3 = min(., -16.4709 + 25.2599*ln(`score2dmpct') + 0.0022*`score2dmpct'^2) if riskreg == 4 & sex == 1
	replace riskage3 = min(., 45.1393 + 12.4337*ln(`score2dmpct') + 0.0823*`score2dmpct'^1) if riskreg == 1 & sex == 2
	replace riskage3 = min(., 44.6316 + 11.3621*ln(`score2dmpct') + 0.0849*`score2dmpct'^1) if riskreg == 2 & sex == 2
	replace riskage3 = min(., 38.9655 + 11.6632*ln(`score2dmpct') + 0.0012*`score2dmpct'^2) if riskreg == 3 & sex == 2
	replace riskage3 = min(., 5.5735 + 18.6020*ln(`score2dmpct') + 0.0017*`score2dmpct'^2) if riskreg == 4 & sex == 2
	replace riskage4 = min(., 46.7859 + 12.1497*ln(`score2dmpct') + 0.0797*`score2dmpct'^1) if riskreg == 1 & sex == 1
	replace riskage4 = min(., 44.9942 + 11.5883*ln(`score2dmpct') + 0.0835*`score2dmpct'^1) if riskreg == 2 & sex == 1
	replace riskage4 = min(., 39.4838 + 12.5391*ln(`score2dmpct') + 0.0935*`score2dmpct'^1) if riskreg == 3 & sex == 1
	replace riskage4 = min(., 2.5341 + 21.1779*ln(`score2dmpct') + 0.0019*`score2dmpct'^2) if riskreg == 4 & sex == 1
	replace riskage4 = min(., 52.8222 + 10.6199*ln(`score2dmpct') + 0.0699*`score2dmpct'^1) if riskreg == 1 & sex == 2
	replace riskage4 = min(., 52.3846 + 9.7108*ln(`score2dmpct') + 0.0720*`score2dmpct'^1) if riskreg == 2 & sex == 2
	replace riskage4 = min(., 47.5868 + 9.9254*ln(`score2dmpct') + 0.0010*`score2dmpct'^2) if riskreg == 3 & sex == 2
	replace riskage4 = min(., 19.2317 + 15.7942*ln(`score2dmpct') + 0.0015*`score2dmpct'^2) if riskreg == 4 & sex == 2
	replace riskage5 = min(., 44.3150 + 12.3421*ln(`score2dmpct') + 0.0838*`score2dmpct'^1) if riskreg == 1 & sex == 1
	replace riskage5 = min(., 42.4866 + 11.7646*ln(`score2dmpct') + 0.0887*`score2dmpct'^1) if riskreg == 2 & sex == 1
	replace riskage5 = min(., 37.2353 + 12.9862*ln(`score2dmpct') + 0.0013*`score2dmpct'^2) if riskreg == 3 & sex == 1
	replace riskage5 = min(., -0.7265 + 21.5636*ln(`score2dmpct') + 0.0019*`score2dmpct'^2) if riskreg == 4 & sex == 1
	replace riskage5 = min(., 50.2173 + 11.2004*ln(`score2dmpct') + 0.0741*`score2dmpct'^1) if riskreg == 1 & sex == 2
	replace riskage5 = min(., 49.7558 + 10.2392*ln(`score2dmpct') + 0.0764*`score2dmpct'^1) if riskreg == 2 & sex == 2
	replace riskage5 = min(., 44.6856 + 10.4804*ln(`score2dmpct') + 0.0011*`score2dmpct'^2) if riskreg == 3 & sex == 2
	replace riskage5 = min(., 14.7232 + 16.6904*ln(`score2dmpct') + 0.0016*`score2dmpct'^2) if riskreg == 4 & sex == 2
	replace riskage6 = min(., 41.7590 + 12.5378*ln(`score2dmpct') + 0.0884*`score2dmpct'^1) if riskreg == 1 & sex == 1
	replace riskage6 = min(., 40.2605 + 12.1714*ln(`score2dmpct') + 0.0012*`score2dmpct'^2) if riskreg == 2 & sex == 1
	replace riskage6 = min(., 34.5610 + 13.2191*ln(`score2dmpct') + 0.0013*`score2dmpct'^2) if riskreg == 3 & sex == 1
	replace riskage6 = min(., -4.1007 + 21.9614*ln(`score2dmpct') + 0.0019*`score2dmpct'^2) if riskreg == 4 & sex == 1
	replace riskage6 = min(., 47.3101 + 11.8478*ln(`score2dmpct') + 0.0788*`score2dmpct'^1) if riskreg == 1 & sex == 2
	replace riskage6 = min(., 46.8223 + 10.8279*ln(`score2dmpct') + 0.0814*`score2dmpct'^1) if riskreg == 2 & sex == 2
	replace riskage6 = min(., 41.4434 + 11.1014*ln(`score2dmpct') + 0.0011*`score2dmpct'^2) if riskreg == 3 & sex == 2
	replace riskage6 = min(., 9.6800 + 17.6941*ln(`score2dmpct') + 0.0016*`score2dmpct'^2) if riskreg == 4 & sex == 2
	foreach i of numlist 1(1)6 {
		replace riskage`i' = 0 if riskage`i' < 0
	}
end


exit		/* force exit to ignore example code below */

* examples
* CVD incidence risk calculations
score2dmrisk
score2dmrisk, riskreg(1) replace
score2dmrisk, percent
score2dmrisk, riskreg(1) replace percent

* using real data
* NB: unavailable variables are generated for illustrative purposes only
webuse nhanes2, clear
gen ccode = "USA"
clonevar ages = age
clonevar hxdiabbin = diabetes
generate smallbin = rbinomial(1, 0.20)
clonevar sbp = bpsystol
generate tchol = tcresult*0.02586
generate hdl = hdresult*0.02586
generate agediab = rnormal(58, 12) if hxdiabbin == 1
generate hba1c_ifcc = rnormal(58, 17) if hxdiabbin == 1
generate egfr_ckdepi = rnormal(74, 20) if hxdiabbin == 1
generate lnegfr_ckdepi = ln(egfr_ckdepi)
score2dmrisk
score2dmrisk, riskreg(1) replace
score2dmrisk, riskreg(2) replace

* managing generated variables
des _SCORE2_START - _SCORE2_END
drop _SCORE2_START - _SCORE2_END
