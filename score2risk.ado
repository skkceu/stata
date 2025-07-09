*! Version : 0.00 12Jul2020
*! Authors : Dr Stephen Kaptoge
*!				 Dr Lisa Pennells
*!				 Steven Hageman
*! Address : Cardiovascular Epidemiology Unit
*!           Dept of Public Health and Primary Care
*!           University of Cambridge
*!           Strangeways Research Laboratory
*!           Worts Causeway
*!           Cambridge CB1 8RN
*!           UK

/*----------------------------------------------------------------------------------------------------------
S Kaptoge: Jul 2020
* Stata program for calculation of SCORE2 CVD risk score, with country/region specific recalibration.

ERFC-derived SCORE2 risk prediction model for estimation of 10-year CVD risk,
with adjustment for non-CVD death as competing risk using IPW weighted 
stratified Cox regression, and reclaibration to four ESC risk regions 

* required variables: 
tchol = total cholesterol in mmol/l
sbp = systolic blood pressure in mmHg
ages = baseline age in years
sex = sex (1 = male, 2 = female)
smallbin = smoking as a binary variable (current = 1, other = 0)
hxdiabbin = history of diabetes as a binary variable (definite = 1, other = 0)

riskreg =  the risk region (1 = low, 2 = Moderate, 3 = high, 4 = very high) 
that applies to the country of residence of the individuals in the dataset 
(see separate table for risk region classification of each country)

relevant descriptive notes.
Jul 2020:	Inital script coded as stata dofile for risk score calculation in validation data
May 2021:	Above converted to current wrapper program for general use
Sep 2021:	Added optional calculation of risk age
Oct 2021:	Added optional calculation of SCORE2-OP risk and merged version SCORE2-OV
Mar 2022:	Added option multdiab(integer 1) to optionally include/exclude consideration of diabetes status
Apr 2022:	Added option nobinary to treat smoking and diabetes as continuous variables if using prevalence
Apr 2022:	Added optional calculation of SCORE2-DM risk (work in progress)
Apr 2022:	Added [if] option to be applied at the end of calculations
Jan 2023:	Added Asia recalibration (work in progress)
------------------------------------------------------------------------------------------------------------*/

capture program drop score2risk
program define score2risk, rclass
	
	version 13
	syntax [if] [, CCODEs riskreg(numlist integer max=1 >=1 <=104) replace percent shortlab rounded(integer 0) ///
	dppscale(real 0) noWEB WEBVERSion(integer 2) DETAILs RISKAGE SCORE2op SCORE2dm SCORE1 MULTDIAB(integer 1) ///
	MULTEXTRA(passthru) FULLADJ(passthru) nobinary noEXTRA]
	
	marksample touse, novarlist
	
	* preserve data in case calculation fails
	preserve
	* display country code mappings to SCORE2 risk regions
	if "`ccodes'" == "ccodes" {
		score2regs, `web' webversion(`webversion') `extra'
		local expvars = "sex ages sbp tchol hdl smallbin hxdiabbin"
		capture confirm numeric variable `expvars'
		local extra = cond(_rc == 0, "`extra'", "noextra")
		if "`extra'" == "noextra" {
			noi di as text "Done for now ..."
			exit
		}
	}
	* web settings
	local webopt = cond("`web'" == "noweb", 0, `webversion')
	local quinoi = "qui"
	`quinoi' {		/* start of quinoi */
		* drop variables from previous calculations if they exist
		capture ds _SCORE2_START - _SCORE2_END
		if _rc == 0 {
			drop `r(varlist)'
		}
		
		* check if required predictor variables exist
		local expvars = "sex ages sbp tchol hdl smallbin hxdiabbin"
		noi confirm numeric variable `expvars'
		capture confirm variable riskreg, exact
		if _rc~=0 {
			if missing("`riskreg'") {
				noi di _newline as text "Please enter SCORE2 risk region (1 = Low, 2 = Moderate, 3 = High, 4 = Very high)"
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
			local orderreg = "order riskreg, after(_SCORE2_START)"
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
				
		* generate marker for start of SCORE2 added variables
		gen _SCORE2_START = .
		`orderreg'		/* order risk region variable if generated above */
		
		* drop existing derived variables
		local dropvars1 = "cages csbp cbmi ctchol chdl"
		local dropvars2 = "lp_m1_ep_cv s010_m1_ep_cv scale1 scale2 cal0_m1_risk cal1_m1_risk"
		local dropvars3 = "score2 score2ov"
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
		local expvars = "ages sbp tchol hdl"
		foreach expvar of varlist `expvars' {
			local expvarlab: variable label `expvar'
			local cvar = "c`expvar'"
			local cvareqn = "`eqn_`cvar''"
			gen `cvar' = `cvareqn'
			label variable `cvar' "Centered `expvarlab': `cvareqn'"
		}
		local hxdiabbin = cond(inlist("`binary'", "nobinary"), "c.hxdiabbin", "1.hxdiabbin")
		local smallbin = cond(inlist("`binary'", "nobinary"), "c.smallbin", "1.smallbin")
		
		* generate linear predictor and baseline survival from CVD Cox FGR, 
		*sex = 1 (Men)
		generat lp_m1_ep_cv =  + 0.3742*(cages) + 0.6012*(`smallbin') + 0.2777*(csbp) + 0.6457*(`hxdiabbin')*`multdiab' + 0.1458*(ctchol) + -0.2698*(chdl) + -0.0755*(`smallbin'#c.cages) + -0.0255*(c.cages#c.csbp) + -0.0983*(`hxdiabbin'#c.cages)*`multdiab' + -0.0281*(c.cages#c.ctchol) + 0.0426*(c.cages#c.chdl) if sex == 1
		generat s010_m1_ep_cv = 0.9605 if sex == 1
		
		*sex = 2 (Women)
		replace lp_m1_ep_cv =  + 0.4648*(cages) + 0.7744*(`smallbin') + 0.3131*(csbp) + 0.8096*(`hxdiabbin')*`multdiab' + 0.1002*(ctchol) + -0.2606*(chdl) + -0.1088*(`smallbin'#c.cages) + -0.0277*(c.cages#c.csbp) + -0.1272*(`hxdiabbin'#c.cages)*`multdiab' + -0.0226*(c.cages#c.ctchol) + 0.0613*(c.cages#c.chdl) if sex == 2
		replace s010_m1_ep_cv = 0.9776 if sex == 2
		
		* replace with rounded coefficients as published for sensitivity comparisons 
		if `rounded' == 1 {
			*sex = 1 (Men)
			replace lp_m1_ep_cv =  ln(1.45)*(cages) + ln(1.82)*(`smallbin') + ln(1.32)*(csbp) + ln(1.91)*(`hxdiabbin')*`multdiab' + ln(1.16)*(ctchol) + ln(0.76)*(chdl) + ln(0.93)*(`smallbin'#c.cages) + ln(0.98)*(c.cages#c.csbp) + ln(0.91)*(`hxdiabbin'#c.cages)*`multdiab' + ln(0.97)*(c.cages#c.ctchol) + ln(1.04)*(c.cages#c.chdl) if sex == 1
			replace s010_m1_ep_cv = 0.9605 if sex == 1
			
			*sex = 2 (Women)
			replace lp_m1_ep_cv =  ln(1.59)*(cages) + ln(2.17)*(`smallbin') + ln(1.37)*(csbp) + ln(2.25)*(`hxdiabbin')*`multdiab' + ln(1.11)*(ctchol) + ln(0.77)*(chdl) + ln(0.90)*(`smallbin'#c.cages) + ln(0.97)*(c.cages#c.csbp) + ln(0.88)*(`hxdiabbin'#c.cages)*`multdiab' + ln(0.98)*(c.cages#c.ctchol) + ln(1.06)*(c.cages#c.chdl) if sex == 2
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
		
		**** Asia: Work in progress ****
		* Asia low risk region
		replace scale1 = -0.323568 if riskreg==101 & sex==1
		replace scale2 =  0.673304 if riskreg==101 & sex==1
		replace scale1 = -0.573924 if riskreg==101 & sex==2
		replace scale2 =  0.601509 if riskreg==101 & sex==2
		
		* Asia moderate risk region
		replace scale1 = -0.323568 if riskreg==102 & sex==1
		replace scale2 =  0.673304 if riskreg==102 & sex==1
		replace scale1 = -0.573924 if riskreg==102 & sex==2
		replace scale2 =  0.601509 if riskreg==102 & sex==2
		
		* Asia high risk region
		replace scale1 = 1.652122 if riskreg==103 & sex==1
		replace scale2 = 0.933416 if riskreg==103 & sex==1
		replace scale1 = 1.561508 if riskreg==103 & sex==2
		replace scale2 = 0.802863 if riskreg==103 & sex==2
		
		* Asia very high risk region
		replace scale1 = 2.109450 if riskreg==104 & sex==1
		replace scale2 = 0.866315 if riskreg==104 & sex==1
		replace scale1 = 2.056977 if riskreg==104 & sex==2
		replace scale2 = 0.735817 if riskreg==104 & sex==2
		
		* replace with rounded coefficients as published for sensitivity comparisons 
		if `rounded' == 1 | `dppscale' > 0 {
			local dpprec = cond(`rounded' == 1 & `dppscale' == 0, 0.01, `dppscale')
			replace scale1 = round(scale1, `dpprec')
			replace scale2 = round(scale2, `dpprec')
		}
		
		* recalibrate the risk estimates for all individuals:
		gen cal1_m1_risk = 1-exp(-exp(scale1 + scale2*ln(-ln(1 - cal0_m1_risk))))
		clonevar score2 = cal1_m1_risk
		
		* label variables informatively
		label variable lp_m1_ep_cv "Linear predictor CVD"
		label variable s010_m1_ep_cv "ERFC baseline survior at 10 years"
		label variable cal0_m1_risk "ERFC 10-yr CVD risk"
		label variable scale1 "SCORE2 rescaling factor 1 (intercept)"
		label variable scale2 "SCORE2 rescaling factor 2 (slope)"
		label variable cal1_m1_risk "SCORE2 10-yr CVD risk"
		label variable score2 "SCORE2 10-yr CVD risk"
		if "`shortlab'" == "shortlab" {
			label variable lp_m1_ep_cv "LP-CVD"
			label variable s010_m1_ep_cv "ERFC-s010"
			label variable cal0_m1_risk "ERFC-CVD"
			label variable scale1 "SCORE2-b0"
			label variable scale2 "SCORE2-b1"
			label variable cal1_m1_risk "SCORE2-CVD"
			label variable score2 "SCORE2-CVD"
		}
		
		if "`percent'" == "percent" {
			local riskvars = "cal?_m?_risk score2"
			foreach var of varlist `riskvars' {
				local varlab: variable label `var'
				replace `var' = 100*`var'		/* convert to percent risk */
				label variable `var' "`varlab' (%)"
			}
		}
		
		* calculate risk age if requested
		if !missing("`riskage'") {
			score2riskage if !missing(score2)		/* call wrapper program below to calculate risk age */
		}

		* generate marker for end of SCORE2 added variables
		gen _SCORE2_END = .
	}		/* end of quinoi */
	
	* cancel restore if successful
	restore, not
	
	* keep minimal variables in dataset if nodetails
	if inlist(`webopt', 1, 2) {
		qui ds _SCORE2_START - _SCORE2_END
		local dropvars = r(varlist)
		local calvars = cond(inlist(`webopt', 1), "cal?_m?_risk score2", "score2")
		if !missing("`riskage'") {
			local calvars = "`calvars' ages riskage?"
		}
		if "`details'" ~= "details" {
			local keepvars = "_SCORE2_START riskreg `calvars' _SCORE2_END"
			qui ds `keepvars'
			local keepvars = r(varlist)
			local dropvars: list dropvars - keepvars
			drop `dropvars'
		}
	}
	
	* calculate SCORE2-OP, SCORE2-DM, or SCORE1 if requested
	if !missing("`score2op'`score2dm'`score1'") {
		qui ds _SCORE2_START - _SCORE2_END
		local score2vars = "`r(varlist)'"
		local exclvars = "_SCORE2_START score2 _SCORE2_END"
		local score2vars: list score2vars - exclvars
		if !missing("`score2vars'") {
			rename (`score2vars') score2_=
		}
		if !missing("`score2op'") {
			qui score2oprisk, `ccodes' riskreg(`riskreg') `replace' `percent' `shortlab' rounded(`rounded') ///
			dppscale(`dppscale') `web' webversion(`webversion') `details' `riskage' multdiab(`multdiab') `binary'
			qui ds _SCORE2OP_START - _SCORE2OP_END
			local score2opvars = "`r(varlist)'"
			local exclvars = "_SCORE2OP_START score2op _SCORE2OP_END"
			local score2opvars: list score2opvars - exclvars
			if !missing("`score2opvars'") {
				rename (`score2opvars') score2op_=
			}
			local sfixes = "lp_m1_ep_cv s010_m1_ep_cv cal0_m1_risk scale1 scale2 cal1_m1_risk"
			foreach sfix in `sfixes' {
				local score2var = "score2_`sfix'"
				local score2opvar = "score2op_`sfix'"
				local score2ovvar = "score2ov_`sfix'"
				capture confirm variable `score2opvar'
				if _rc == 0 {
					gen `score2ovvar' = cond(ages < 70, `score2var', `score2opvar')
					local varlab: variable label `score2opvar'
					local varlab: subinstr local varlab "CONOR" "ERFC/CONOR"
					local varlab: subinstr local varlab "SCORE2-OP" "SCORE2-OV"
					label variable `score2ovvar' "`varlab'"
				}
			}
			gen score2ov = cond(ages < 70, score2, score2op)
			local varlab: variable label score2op
			local varlab: subinstr local varlab "SCORE2-OP" "SCORE2-OV"
			label variable score2ov "`varlab'"
		}
		if !missing("`score2dm'") {
			qui score2dmrisk, `ccodes' riskreg(`riskreg') `replace' `percent' `shortlab' rounded(`rounded') ///
			dppscale(`dppscale') `web' webversion(`webversion') `details' `riskage' multdiab(`multdiab') ///
			`multextra' `fulladj' `binary'
			qui ds _SCORE2DM_START - _SCORE2DM_END
			local score2dmvars = "`r(varlist)'"
			local exclvars = "_SCORE2DM_START score2dm _SCORE2DM_END"
			local score2dmvars: list score2dmvars - exclvars
			if !missing("`score2dmvars'") {
				rename (`score2dmvars') score2dm_=
			}
		}
		if !missing("`score1'") {
			qui score1risk, `ccodes' score1rreg(`riskreg') `replace' `percent' `shortlab' rounded(`rounded') ///
			dppscale(`dppscale') `web' webversion(`webversion') `details' `riskage' multdiab(`multdiab')
			qui ds _SCORE1_START - _SCORE1_END
			local score1vars = "`r(varlist)'"
			local exclvars = "_SCORE1_START score1 _SCORE1_END"
			local score1vars: list score1vars - exclvars
			if !missing("`score1vars'") {
				rename (`score1vars') score1_=
			}
		}
		* calculate risk age if requested
		if !missing("`score2op'") & !missing("`riskage'") {
			qui score2ovriskage if !missing(score2ov)		/* call wrapper program below to calculate risk age */
			local score2ovvars = "riskage?"
		}
		if !missing("`score2ovvars'") {
			rename (`score2ovvars') score2ov_=
			local score2ovvars = "score2ov_*"
		}
		if !missing("`score2op'") {
			local score2ovvars = "score2ov* score2ov `score2ovvars'"
			order `score2ovvars', before(_SCORE2OP_END)
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
	qui recode _SCORE2_START - _SCORE2_END (* = .) if !(`touse')
	des _SCORE2_START - _SCORE2_END
	if !missing("`score2op'") {
		qui recode _SCORE2OP_START - _SCORE2OP_END (* = .) if !(`touse')
		des _SCORE2OP_START - _SCORE2OP_END
	}
	if !missing("`score2dm'") {
		qui recode _SCORE2DM_START - _SCORE2DM_END (* = .) if !(`touse')
		des _SCORE2DM_START - _SCORE2DM_END
	}
	if !missing("`score1'") {
		qui recode _SCORE1_START - _SCORE1_END (* = .) if !(`touse')
		des _SCORE1_START - _SCORE1_END
	}
	preserve
	qui keep if !missing(score2)
	local calvars = cond(inlist(`webopt', 0, 1), "cal?_m?_risk score2", "score2")
	if !missing("`score2op'`score2dm'") {
		local calvars = cond(inlist(`webopt', 0, 1), "score2*_cal?_m?_risk score2 score2??", "score2 score2??")
	}
	if !missing("`score1'") {
		local calvars = "`calvars' score1"
	}
	if !missing("`riskage'") {
		local calvars = "`calvars' ages *riskage?"
	}
	capture confirm variable riskreg, exact
	if _rc ~= 0 {
		clonevar riskreg = score2_riskreg
	}
	bysort riskreg: tabstat `calvars', by(sex) stats(n mean sd min p25 p50 p75 max) ///
	col(stats) longstub varw(16)
	if !missing("`riskage'") {
		notes `calvars'
	}
	restore
	
end

* wrapper program to calculate risk age for SCORE2
capture program drop score2riskage
program define score2riskage, rclass
	
	version 13
	syntax [if] [in] [, *]
	marksample touse, novarlist
	
	local dropvars = "riskage1 riskage2 riskage3 riskage4 riskage5 riskage6"
	foreach var in `dropvars' {
		capture drop `var'
	}
	foreach i of numlist 1(1)6 {
		gen riskage`i' = .
		label variable riskage`i' "SCORE2 risk age`i'"
	}
	notes riskage1: "SCORE2 risk age1: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage2: "SCORE2 risk age2: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage3: "SCORE2 risk age3: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage4: "SCORE2 risk age4: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage5: "SCORE2 risk age5: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage6: "SCORE2 risk age6: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	tempvar score2pct
	summ score2
	gen `score2pct' = cond(r(max) < 1, score2*100, score2)
	* risk age estimation equations derived separately (see score2risk_test.do)
	replace riskage1 = min(., 31.8020 + 17.0321*ln(`score2pct') + 0.1062*`score2pct'^1) if riskreg == 1 & sex == 1
	replace riskage1 = min(., 29.7734 + 15.8604*ln(`score2pct') + 0.1085*`score2pct'^1) if riskreg == 2 & sex == 1
	replace riskage1 = min(., 33.7703 + 13.8287*ln(`score2pct') + 0.0014*`score2pct'^2) if riskreg == 3 & sex == 1
	replace riskage1 = min(., 20.9922 + 15.6804*ln(`score2pct') + 0.0015*`score2pct'^2) if riskreg == 4 & sex == 1
	replace riskage1 = min(., 42.2955 + 14.4114*ln(`score2pct') + 0.0885*`score2pct'^1) if riskreg == 1 & sex == 2
	replace riskage1 = min(., 41.6400 + 13.1128*ln(`score2pct') + 0.0882*`score2pct'^1) if riskreg == 2 & sex == 2
	replace riskage1 = min(., 42.4021 + 10.9250*ln(`score2pct') + 0.0011*`score2pct'^2) if riskreg == 3 & sex == 2
	replace riskage1 = min(., 30.8637 + 12.3014*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 4 & sex == 2
	replace riskage2 = min(., 33.1255 + 16.7621*ln(`score2pct') + 0.1040*`score2pct'^1) if riskreg == 1 & sex == 1
	replace riskage2 = min(., 31.1281 + 15.6116*ln(`score2pct') + 0.1061*`score2pct'^1) if riskreg == 2 & sex == 1
	replace riskage2 = min(., 35.0711 + 13.6005*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 3 & sex == 1
	replace riskage2 = min(., 22.5057 + 15.4209*ln(`score2pct') + 0.0015*`score2pct'^2) if riskreg == 4 & sex == 1
	replace riskage2 = min(., 43.3235 + 14.2135*ln(`score2pct') + 0.0869*`score2pct'^1) if riskreg == 1 & sex == 2
	replace riskage2 = min(., 42.6775 + 12.9342*ln(`score2pct') + 0.0864*`score2pct'^1) if riskreg == 2 & sex == 2
	replace riskage2 = min(., 43.4322 + 10.7710*ln(`score2pct') + 0.0011*`score2pct'^2) if riskreg == 3 & sex == 2
	replace riskage2 = min(., 32.0562 + 12.1299*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 4 & sex == 2
	replace riskage3 = min(., 34.4073 + 16.5005*ln(`score2pct') + 0.1020*`score2pct'^1) if riskreg == 1 & sex == 1
	replace riskage3 = min(., 32.4405 + 15.3703*ln(`score2pct') + 0.1038*`score2pct'^1) if riskreg == 2 & sex == 1
	replace riskage3 = min(., 35.9987 + 13.1132*ln(`score2pct') + 0.1002*`score2pct'^1) if riskreg == 3 & sex == 1
	replace riskage3 = min(., 23.9706 + 15.1698*ln(`score2pct') + 0.0014*`score2pct'^2) if riskreg == 4 & sex == 1
	replace riskage3 = min(., 44.3235 + 14.0210*ln(`score2pct') + 0.0854*`score2pct'^1) if riskreg == 1 & sex == 2
	replace riskage3 = min(., 43.6867 + 12.7603*ln(`score2pct') + 0.0847*`score2pct'^1) if riskreg == 2 & sex == 2
	replace riskage3 = min(., 44.4338 + 10.6213*ln(`score2pct') + 0.0011*`score2pct'^2) if riskreg == 3 & sex == 2
	replace riskage3 = min(., 33.2161 + 11.9629*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 4 & sex == 2
	replace riskage4 = min(., 29.7912 + 17.6645*ln(`score2pct') + 0.1099*`score2pct'^1) if riskreg == 1 & sex == 1
	replace riskage4 = min(., 27.6912 + 16.4483*ln(`score2pct') + 0.1121*`score2pct'^1) if riskreg == 2 & sex == 1
	replace riskage4 = min(., 31.8112 + 14.3531*ln(`score2pct') + 0.0014*`score2pct'^2) if riskreg == 3 & sex == 1
	replace riskage4 = min(., 18.5437 + 16.2806*ln(`score2pct') + 0.0015*`score2pct'^2) if riskreg == 4 & sex == 1
	replace riskage4 = min(., 41.3630 + 14.7492*ln(`score2pct') + 0.0902*`score2pct'^1) if riskreg == 1 & sex == 2
	replace riskage4 = min(., 40.6942 + 13.4203*ln(`score2pct') + 0.0897*`score2pct'^1) if riskreg == 2 & sex == 2
	replace riskage4 = min(., 41.4664 + 11.1860*ln(`score2pct') + 0.0011*`score2pct'^2) if riskreg == 3 & sex == 2
	replace riskage4 = min(., 29.6553 + 12.5985*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 4 & sex == 2
	replace riskage5 = min(., 31.1962 + 17.3744*ln(`score2pct') + 0.1076*`score2pct'^1) if riskreg == 1 & sex == 1
	replace riskage5 = min(., 29.1295 + 16.1809*ln(`score2pct') + 0.1095*`score2pct'^1) if riskreg == 2 & sex == 1
	replace riskage5 = min(., 32.8820 + 13.8010*ln(`score2pct') + 0.1056*`score2pct'^1) if riskreg == 3 & sex == 1
	replace riskage5 = min(., 20.1551 + 16.0008*ln(`score2pct') + 0.0015*`score2pct'^2) if riskreg == 4 & sex == 1
	replace riskage5 = min(., 42.4279 + 14.5420*ln(`score2pct') + 0.0886*`score2pct'^1) if riskreg == 1 & sex == 2
	replace riskage5 = min(., 41.7688 + 13.2333*ln(`score2pct') + 0.0879*`score2pct'^1) if riskreg == 2 & sex == 2
	replace riskage5 = min(., 42.5337 + 11.0244*ln(`score2pct') + 0.0011*`score2pct'^2) if riskreg == 3 & sex == 2
	replace riskage5 = min(., 30.8930 + 12.4184*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 4 & sex == 2
	replace riskage6 = min(., 32.5555 + 17.0934*ln(`score2pct') + 0.1054*`score2pct'^1) if riskreg == 1 & sex == 1
	replace riskage6 = min(., 30.5213 + 15.9218*ln(`score2pct') + 0.1071*`score2pct'^1) if riskreg == 2 & sex == 1
	replace riskage6 = min(., 34.2126 + 13.5833*ln(`score2pct') + 0.1030*`score2pct'^1) if riskreg == 3 & sex == 1
	replace riskage6 = min(., 21.7128 + 15.7304*ln(`score2pct') + 0.0015*`score2pct'^2) if riskreg == 4 & sex == 1
	replace riskage6 = min(., 43.4630 + 14.3405*ln(`score2pct') + 0.0870*`score2pct'^1) if riskreg == 1 & sex == 2
	replace riskage6 = min(., 42.8136 + 13.0512*ln(`score2pct') + 0.0862*`score2pct'^1) if riskreg == 2 & sex == 2
	replace riskage6 = min(., 43.5709 + 10.8674*ln(`score2pct') + 0.0011*`score2pct'^2) if riskreg == 3 & sex == 2
	replace riskage6 = min(., 32.0959 + 12.2432*ln(`score2pct') + 0.0013*`score2pct'^2) if riskreg == 4 & sex == 2
	foreach i of numlist 1(1)6 {
		replace riskage`i' = 0 if riskage`i' < 0
	}
end

* wrapper program to calculate risk age for SCORE2-OV
capture program drop score2ovriskage
program define score2ovriskage, rclass
	
	version 13
	syntax [if] [in] [, *]
	marksample touse, novarlist
	
	local dropvars = "riskage1 riskage2 riskage3 riskage4 riskage5 riskage6"
	foreach var in `dropvars' {
		capture drop `var'
	}
	foreach i of numlist 1(1)6 {
		gen riskage`i' = .
		label variable riskage`i' "SCORE2-OV risk age`i'"
	}
	notes riskage1: "SCORE2-OV risk age1: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage2: "SCORE2-OV risk age2: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage3: "SCORE2-OV risk age3: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage4: "SCORE2-OV risk age4: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage5: "SCORE2-OV risk age5: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage6: "SCORE2-OV risk age6: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	tempvar score2ovpct
	summ score2ov
	gen `score2ovpct' = cond(r(max) < 1, score2ov*100, score2ov)
	* risk age estimation equations derived separately (see score2ovrisk_test.do)
	replace riskage1 = min(., 32.1149 + 17.3548*ln(`score2ovpct') + -0.0270*`score2ovpct'^1) if riskreg == 1 & sex == 1
	replace riskage1 = min(., 29.8873 + 16.1309*ln(`score2ovpct') + 0.1140*ln(`score2ovpct')) if riskreg == 2 & sex == 1
	replace riskage1 = min(., 33.3984 + 13.5090*ln(`score2ovpct') + 0.1361*`score2ovpct'^1) if riskreg == 3 & sex == 1
	replace riskage1 = min(., 21.0132 + 15.5696*ln(`score2ovpct') + 0.0001*`score2ovpct'^3) if riskreg == 4 & sex == 1
	replace riskage1 = min(., 42.6782 + 14.6260*ln(`score2ovpct') + -0.0524*`score2ovpct'^1) if riskreg == 1 & sex == 2
	replace riskage1 = min(., 41.9906 + 13.3091*ln(`score2ovpct') + -0.0065*`score2ovpct'^1) if riskreg == 2 & sex == 2
	replace riskage1 = min(., 42.0921 + 10.7703*ln(`score2ovpct') + 0.0925*`score2ovpct'^1) if riskreg == 3 & sex == 2
	replace riskage1 = min(., 30.7619 + 12.1550*ln(`score2ovpct') + 0.0021*`score2ovpct'^2) if riskreg == 4 & sex == 2
	replace riskage2 = min(., 33.4867 + 17.1134*ln(`score2ovpct') + -0.0401*`score2ovpct'^1) if riskreg == 1 & sex == 1
	replace riskage2 = min(., 30.2634 + 0.9744*`score2ovpct'^-.5 + 16.4421*ln(`score2ovpct')) if riskreg == 2 & sex == 1
	replace riskage2 = min(., 34.7308 + 13.3293*ln(`score2ovpct') + 0.1240*`score2ovpct'^1) if riskreg == 3 & sex == 1
	replace riskage2 = min(., 22.5159 + 15.3055*ln(`score2ovpct') + 0.0001*`score2ovpct'^3) if riskreg == 4 & sex == 1
	replace riskage2 = min(., 43.7244 + 14.4280*ln(`score2ovpct') + -0.0588*`score2ovpct'^1) if riskreg == 1 & sex == 2
	replace riskage2 = min(., 43.0469 + 13.1303*ln(`score2ovpct') + -0.0120*`score2ovpct'^1) if riskreg == 2 & sex == 2
	replace riskage2 = min(., 43.1369 + 10.6292*ln(`score2ovpct') + 0.0869*`score2ovpct'^1) if riskreg == 3 & sex == 2
	replace riskage2 = min(., 31.9370 + 11.9763*ln(`score2ovpct') + 0.0020*`score2ovpct'^2) if riskreg == 4 & sex == 2
	replace riskage3 = min(., 34.8172 + 16.8751*ln(`score2ovpct') + -0.0523*`score2ovpct'^1) if riskreg == 1 & sex == 1
	replace riskage3 = min(., 32.4457 + 0.1623*`score2ovpct'^-1 + 15.9291*ln(`score2ovpct')) if riskreg == 2 & sex == 1
	replace riskage3 = min(., 36.0222 + 13.1506*ln(`score2ovpct') + 0.1126*`score2ovpct'^1) if riskreg == 3 & sex == 1
	replace riskage3 = min(., 23.9692 + 15.0476*ln(`score2ovpct') + 0.0000*`score2ovpct'^3) if riskreg == 4 & sex == 1
	replace riskage3 = min(., 44.7411 + 14.2343*ln(`score2ovpct') + -0.0649*`score2ovpct'^1) if riskreg == 1 & sex == 2
	replace riskage3 = min(., 44.0737 + 12.9553*ln(`score2ovpct') + -0.0172*`score2ovpct'^1) if riskreg == 2 & sex == 2
	replace riskage3 = min(., 44.5389 + 10.6740*ln(`score2ovpct') + 0.0009*`score2ovpct'^2) if riskreg == 3 & sex == 2
	replace riskage3 = min(., 33.0772 + 11.8015*ln(`score2ovpct') + 0.0020*`score2ovpct'^2) if riskreg == 4 & sex == 2
	replace riskage4 = min(., 30.1749 + 18.0442*ln(`score2ovpct') + -0.0832*`score2ovpct'^1) if riskreg == 1 & sex == 1
	replace riskage4 = min(., 27.9859 + 16.8316*ln(`score2ovpct') + -0.0143*`score2ovpct'^1) if riskreg == 2 & sex == 1
	replace riskage4 = min(., 31.5567 + 14.0803*ln(`score2ovpct') + 0.0960*`score2ovpct'^1) if riskreg == 3 & sex == 1
	replace riskage4 = min(., 18.5819 + 16.1827*ln(`score2ovpct') + 0.0000*`score2ovpct'^3) if riskreg == 4 & sex == 1
	replace riskage4 = min(., 41.7339 + 14.9573*ln(`score2ovpct') + -0.0644*`score2ovpct'^1) if riskreg == 1 & sex == 2
	replace riskage4 = min(., 41.2309 + 13.6860*ln(`score2ovpct') + -0.1763*`score2ovpct'^.5) if riskreg == 2 & sex == 2
	replace riskage4 = min(., 41.5470 + 11.2354*ln(`score2ovpct') + 0.0010*`score2ovpct'^2) if riskreg == 3 & sex == 2
	replace riskage4 = min(., 29.5640 + 12.4382*ln(`score2ovpct') + 0.0021*`score2ovpct'^2) if riskreg == 4 & sex == 2
	replace riskage5 = min(., 31.6146 + 17.7672*ln(`score2ovpct') + -0.0926*`score2ovpct'^1) if riskreg == 1 & sex == 1
	replace riskage5 = min(., 29.4631 + 16.5781*ln(`score2ovpct') + -0.0231*`score2ovpct'^1) if riskreg == 2 & sex == 1
	replace riskage5 = min(., 32.9580 + 13.8741*ln(`score2ovpct') + 0.0859*`score2ovpct'^1) if riskreg == 3 & sex == 1
	replace riskage5 = min(., 20.1861 + 15.8901*ln(`score2ovpct') + 0.0000*`score2ovpct'^3) if riskreg == 4 & sex == 1
	replace riskage5 = min(., 42.8140 + 14.7488*ln(`score2ovpct') + -0.0705*`score2ovpct'^1) if riskreg == 1 & sex == 2
	replace riskage5 = min(., 42.3909 + 13.5209*ln(`score2ovpct') + -0.2360*`score2ovpct'^.5) if riskreg == 2 & sex == 2
	replace riskage5 = min(., 42.6066 + 11.0646*ln(`score2ovpct') + 0.0009*`score2ovpct'^2) if riskreg == 3 & sex == 2
	replace riskage5 = min(., 30.7834 + 12.2505*ln(`score2ovpct') + 0.0020*`score2ovpct'^2) if riskreg == 4 & sex == 2
	replace riskage6 = min(., 33.0085 + 17.4959*ln(`score2ovpct') + -0.1013*`score2ovpct'^1) if riskreg == 1 & sex == 1
	replace riskage6 = min(., 30.8942 + 16.3293*ln(`score2ovpct') + -0.0312*`score2ovpct'^1) if riskreg == 2 & sex == 1
	replace riskage6 = min(., 34.3143 + 13.6705*ln(`score2ovpct') + 0.0764*`score2ovpct'^1) if riskreg == 3 & sex == 1
	replace riskage6 = min(., 21.7347 + 15.6055*ln(`score2ovpct') + 0.0000*`score2ovpct'^3) if riskreg == 4 & sex == 1
	replace riskage6 = min(., 43.8630 + 14.5450*ln(`score2ovpct') + -0.0764*`score2ovpct'^1) if riskreg == 1 & sex == 2
	replace riskage6 = min(., 43.5170 + 13.3578*ln(`score2ovpct') + -0.2928*`score2ovpct'^.5) if riskreg == 2 & sex == 2
	replace riskage6 = min(., 43.6337 + 10.8983*ln(`score2ovpct') + 0.0009*`score2ovpct'^2) if riskreg == 3 & sex == 2
	replace riskage6 = min(., 31.9659 + 12.0671*ln(`score2ovpct') + 0.0019*`score2ovpct'^2) if riskreg == 4 & sex == 2
	foreach i of numlist 1(1)6 {
		replace riskage`i' = 0 if riskage`i' < 0
	}
end


exit		/* force exit to ignore example code below */

* examples
* CVD incidence risk calculations
score2risk
score2risk, riskreg(1) replace
score2risk, percent
score2risk, riskreg(1) replace percent

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
score2risk
score2risk, riskreg(1) replace
score2risk, riskreg(2) replace
score2risk, riskreg(2) replace score2op
score2risk, riskreg(2) replace score2dm

* managing generated variables
des _SCORE2_START - _SCORE2_END
drop _SCORE2_START - _SCORE2_END
