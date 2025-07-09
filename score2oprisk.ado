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
* Stata program for calculation of SCORE2-OP CVD risk score, with country/region specific recalibration.

CONOR-derived SCORE2-OP risk prediction model for estimation of 10-year CVD risk,
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
Sep 2021:	Above converted to current wrapper program for general use
Oct 2021:	Added optional calculation of risk age
Mar 2022:	Added option multdiab(integer 1) to optionally include/exclude consideration of diabetes status
Apr 2022:	Added option nobinary to treat smoking and diabetes as continuous variables if using prevalence
Apr 2022:	Added [if] option to be applied at the end of calculations
------------------------------------------------------------------------------------------------------------*/

capture program drop score2oprisk
program define score2oprisk, rclass
	
	version 13
	syntax [if] [, CCODEs riskreg(numlist integer max=1 >=1 <=104) replace percent shortlab rounded(integer 0) ///
	dppscale(real 0) noWEB WEBVERSion(integer 2) DETAILs RISKAGE MULTDIAB(integer 1) nobinary]
	
	marksample touse, novarlist
	
	* preserve data in case calculation fails
	preserve
	* web settings
	local webopt = cond("`web'" == "noweb", 0, `webversion')
	local quinoi = "qui"
	`quinoi' {		/* start of quinoi */
		* drop variables from previous calculations if they exist
		capture ds _SCORE2OP_START - _SCORE2OP_END
		if _rc == 0 {
			drop `r(varlist)'
		}
		
		* check if required predictor variables exist
		local expvars = "sex ages sbp tchol hdl smallbin hxdiabbin"
		noi confirm numeric variable `expvars'
		capture confirm variable riskreg, exact
		if _rc~=0 {
			if missing("`riskreg'") {
				noi di _newline as text "Please enter SCORE2-OP risk region (1 = Low, 2 = Moderate, 3 = High, 4 = Very high)"
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
			local orderreg = "order riskreg, after(_SCORE2OP_START)"
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
				
		* generate marker for start of SCORE2-OP added variables
		gen _SCORE2OP_START = .
		`orderreg'		/* order risk region variable if generated above */
		
		* drop existing derived variables
		local dropvars1 = "cages csbp cbmi ctchol chdl"
		local dropvars2 = "lp_m1_ep_cv s010_m1_ep_cv scale1 scale2 cal0_m1_risk cal1_m1_risk"
		local dropvars3 = "score2op"
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
		local eqn_cages = "(ages - 73)/1"
		local eqn_csbp = "(sbp - 150)/1"
		local eqn_ctchol = "(tchol - 6)/1"
		local eqn_chdl = "(hdl - 1.4)/1"
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
		local lp_mean1 = 0.09285225
		generat lp_m1_ep_cv = 0.0633767915*(cages) + 0.4245248910*(`hxdiabbin')*`multdiab' + -0.0173654217*(`hxdiabbin'#c.cages)*`multdiab' + 0.3523957899*(`smallbin') + ///
									-0.0247205945*(`smallbin'#c.cages) + 0.0094362177*(csbp) + -0.0004794522*(c.cages#c.csbp) + 0.0849594156*(ctchol) + ///
									0.0073336663*(c.cages#c.ctchol) + -0.3563571842*(chdl) + 0.0091317857*(c.cages#c.chdl) - `lp_mean1' if sex == 1
		generat s010_m1_ep_cv = 1 - 0.24235298 if sex == 1
		
		*sex = 2 (Women)
		local lp_mean2 = 0.2287676
		replace lp_m1_ep_cv = 0.0789191229*(cages) + 0.6009874438*(`hxdiabbin')*`multdiab' + -0.0106909701*(`hxdiabbin'#c.cages)*`multdiab' + 0.4921147010*(`smallbin') + ///
									-0.0255496313*(`smallbin'#c.cages) + 0.0101993091*(csbp) + -0.0004404728*(c.cages#c.csbp) + 0.0604939598*(ctchol) + ///
									-0.0009121461*(c.cages#c.ctchol) + -0.3040280248*(chdl) + 0.0153779923*(c.cages#c.chdl) - `lp_mean2' if sex == 2
		replace s010_m1_ep_cv = 1 - 0.19178759 if sex == 2

		* replace with rounded coefficients as published for sensitivity comparisons 
		if `rounded' == 1 {
			*sex = 1 (Men)
			local lp_mean1 = 0.0929
			replace lp_m1_ep_cv = 0.0634*(cages) + 0.4245*(`hxdiabbin')*`multdiab' + -0.0174*(`hxdiabbin'#c.cages)*`multdiab' + 0.3524*(`smallbin') + ///
										-0.0247*(`smallbin'#c.cages) + 0.0094*(csbp) + -0.0005*(c.cages#c.csbp) + 0.0850*(ctchol) + ///
										0.0073*(c.cages#c.ctchol) + -0.3564*(chdl) + 0.0091*(c.cages#c.chdl) - `lp_mean1' if sex == 1
			replace s010_m1_ep_cv = 1 - 0.2424 if sex == 1
			
			*sex = 2 (Women)
			local lp_mean2 = 0.2288
			replace lp_m1_ep_cv = 0.0789*(cages) + 0.6010*(`hxdiabbin')*`multdiab' + -0.0107*(`hxdiabbin'#c.cages)*`multdiab' + 0.4921*(`smallbin') + ///
										-0.0255*(`smallbin'#c.cages) + 0.0102*(csbp) + -0.0004*(c.cages#c.csbp) + 0.0605*(ctchol) + ///
										-0.0009*(c.cages#c.ctchol) + -0.3040*(chdl) + 0.0154*(c.cages#c.chdl) - `lp_mean2' if sex == 2
			replace s010_m1_ep_cv = 1 - 0.1918 if sex == 2
		}
		
		* Calculate uncalibrated 10-year risk calculation from linear predictors and baseline survival
		gen cal0_m1_risk = 1 - s010_m1_ep_cv^exp(lp_m1_ep_cv)
		
		* generate recalibration factors depending on sex and risk region of the individuals
		gen scale1 = .
		gen scale2 = .

		**** Europe: Published work ****
		* Europe low risk region
		replace scale1 = -0.343484 if riskreg == 1 & sex == 1
		replace scale2 = 1.193038 if riskreg == 1 & sex == 1
		replace scale1 = -0.516054 if riskreg == 1 & sex == 2
		replace scale2 = 1.007238 if riskreg == 1 & sex == 2
		
		* Europe moderate risk region
		replace scale1 = 0.005320 if riskreg == 2 & sex == 1
		replace scale2 = 1.248514 if riskreg == 2 & sex == 1
		replace scale1 = -0.101044 if riskreg == 2 & sex == 2
		replace scale2 = 1.099621 if riskreg == 2 & sex == 2
		
		* Europe high risk region
		replace scale1 = 0.084305 if riskreg == 3 & sex == 1
		replace scale2 = 1.151902 if riskreg == 3 & sex == 1
		replace scale1 = 0.375720 if riskreg == 3 & sex == 2
		replace scale2 = 1.089031 if riskreg == 3 & sex == 2
		
		* Europe very high risk region
		replace scale1 = 0.045959 if riskreg == 4 & sex == 1
		replace scale2 = 0.703625 if riskreg == 4 & sex == 1
		replace scale1 = 0.379293 if riskreg == 4 & sex == 2
		replace scale2 = 0.691530 if riskreg == 4 & sex == 2
		
		* replace with rounded coefficients as published for sensitivity comparisons 
		if `rounded' == 1 | `dppscale' > 0 {
			local dpprec = cond(`rounded' == 1 & `dppscale' == 0, 0.01, `dppscale')
			replace scale1 = round(scale1, `dpprec')
			replace scale2 = round(scale2, `dpprec')
		}
		
		* recalibrate the risk estimates for all individuals:
		gen cal1_m1_risk = 1-exp(-exp(scale1 + scale2*ln(-ln(1 - cal0_m1_risk))))
		clonevar score2op = cal1_m1_risk
		
		* label variables informatively
		label variable lp_m1_ep_cv "Linear predictor CVD"
		label variable s010_m1_ep_cv "CONOR baseline survior at 10 years"
		label variable cal0_m1_risk "CONOR 10-yr CVD risk"
		label variable scale1 "SCORE2-OP rescaling factor 1 (intercept)"
		label variable scale2 "SCORE2-OP rescaling factor 2 (slope)"
		label variable cal1_m1_risk "SCORE2-OP 10-yr CVD risk"
		label variable score2op "SCORE2-OP 10-yr CVD risk"
		if "`shortlab'" == "shortlab" {
			label variable lp_m1_ep_cv "LP-CVD"
			label variable s010_m1_ep_cv "UTR-s010"
			label variable cal0_m1_risk "UTR-CVD"
			label variable scale1 "SCORE2-OP-b0"
			label variable scale2 "SCORE2-OP-b1"
			label variable cal1_m1_risk "SCORE2-OP-CVD"
			label variable score2op "SCORE2-OP-CVD"
		}
		
		if "`percent'" == "percent" {
			local riskvars = "cal?_m?_risk score2op"
			foreach var of varlist `riskvars' {
				local varlab: variable label `var'
				replace `var' = 100*`var'		/* convert to percent risk */
				label variable `var' "`varlab' (%)"
			}
		}
		
		* calculate risk age if requested
		if !missing("`riskage'") {
			score2opriskage if !missing(score2op)		/* call wrapper program below to calculate risk age */
		}

		* generate marker for end of SCORE2 added variables
		gen _SCORE2OP_END = .
	}		/* end of quinoi */
	
	* cancel restore if successful
	restore, not
	
	* keep minimal variables in dataset if nodetails
	if inlist(`webopt', 1, 2) {
		qui ds _SCORE2OP_START - _SCORE2OP_END
		local dropvars = r(varlist)
		local calvars = cond(inlist(`webopt', 1), "cal?_m?_risk score2op", "score2op")
		if !missing("`riskage'") {
			local calvars = "`calvars' ages riskage?"
		}
		if "`details'" ~= "details" {
			local keepvars = "_SCORE2OP_START riskreg `calvars' _SCORE2OP_END"
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
	qui recode _SCORE2OP_START - _SCORE2OP_END (* = .) if !(`touse')
	des _SCORE2OP_START - _SCORE2OP_END
	preserve
	qui keep if !missing(score2)
	local calvars = cond(inlist(`webopt', 0, 1), "cal?_m?_risk score2op", "score2op")
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

* wrapper program to calculate risk age for SCORE2-OP
capture program drop score2opriskage
program define score2opriskage, rclass
	
	version 13
	syntax [if] [in] [, *]
	marksample touse, novarlist
	
	local dropvars = "riskage1 riskage2 riskage3 riskage4 riskage5 riskage6"
	foreach var in `dropvars' {
		capture drop `var'
	}
	foreach i of numlist 1(1)6 {
		gen riskage`i' = .
		label variable riskage`i' "SCORE2-OP risk age`i'"
	}
	notes riskage1: "SCORE2-OP risk age1: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage2: "SCORE2-OP risk age2: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage3: "SCORE2-OP risk age3: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage4: "SCORE2-OP risk age4: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage5: "SCORE2-OP risk age5: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage6: "SCORE2-OP risk age6: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	tempvar score2oppct
	summ score2op
	gen `score2oppct' = cond(r(max) < 1, score2op*100, score2op)
	* risk age estimation equations derived separately (see score2oprisk_test.do)
	replace riskage1 = min(., 39.9846 + 13.8987*ln(`score2oppct') + 0.0910*`score2oppct'^1) if riskreg == 1 & sex == 1
	replace riskage1 = min(., 37.9434 + 13.2507*ln(`score2oppct') + 0.0952*`score2oppct'^1) if riskreg == 2 & sex == 1
	replace riskage1 = min(., 31.6487 + 14.3310*ln(`score2oppct') + 0.1070*`score2oppct'^1) if riskreg == 3 & sex == 1
	replace riskage1 = min(., -10.9389 + 24.3541*ln(`score2oppct') + 0.0021*`score2oppct'^2) if riskreg == 4 & sex == 1
	replace riskage1 = min(., 49.8338 + 11.2458*ln(`score2oppct') + 0.0747*`score2oppct'^1) if riskreg == 1 & sex == 2
	replace riskage1 = min(., 49.3691 + 10.2798*ln(`score2oppct') + 0.0772*`score2oppct'^1) if riskreg == 2 & sex == 2
	replace riskage1 = min(., 44.2790 + 10.5245*ln(`score2oppct') + 0.0011*`score2oppct'^2) if riskreg == 3 & sex == 2
	replace riskage1 = min(., 14.1880 + 16.7619*ln(`score2oppct') + 0.0016*`score2oppct'^2) if riskreg == 4 & sex == 2
	replace riskage2 = min(., 42.0019 + 13.3679*ln(`score2oppct') + 0.0876*`score2oppct'^1) if riskreg == 1 & sex == 1
	replace riskage2 = min(., 40.0359 + 12.7463*ln(`score2oppct') + 0.0917*`score2oppct'^1) if riskreg == 2 & sex == 1
	replace riskage2 = min(., 33.9787 + 13.7874*ln(`score2oppct') + 0.1029*`score2oppct'^1) if riskreg == 3 & sex == 1
	replace riskage2 = min(., -6.8898 + 23.3883*ln(`score2oppct') + 0.0021*`score2oppct'^2) if riskreg == 4 & sex == 1
	replace riskage2 = min(., 50.9631 + 10.9725*ln(`score2oppct') + 0.0729*`score2oppct'^1) if riskreg == 1 & sex == 2
	replace riskage2 = min(., 50.5088 + 10.0309*ln(`score2oppct') + 0.0754*`score2oppct'^1) if riskreg == 2 & sex == 2
	replace riskage2 = min(., 45.5486 + 10.2633*ln(`score2oppct') + 0.0010*`score2oppct'^2) if riskreg == 3 & sex == 2
	replace riskage2 = min(., 16.2140 + 16.3402*ln(`score2oppct') + 0.0015*`score2oppct'^2) if riskreg == 4 & sex == 2
	replace riskage3 = min(., 43.8704 + 12.8760*ln(`score2oppct') + 0.0845*`score2oppct'^1) if riskreg == 1 & sex == 1
	replace riskage3 = min(., 41.9744 + 12.2787*ln(`score2oppct') + 0.0885*`score2oppct'^1) if riskreg == 2 & sex == 1
	replace riskage3 = min(., 36.1376 + 13.2833*ln(`score2oppct') + 0.0992*`score2oppct'^1) if riskreg == 3 & sex == 1
	replace riskage3 = min(., -3.1451 + 22.4955*ln(`score2oppct') + 0.0020*`score2oppct'^2) if riskreg == 4 & sex == 1
	replace riskage3 = min(., 52.0387 + 10.7122*ln(`score2oppct') + 0.0712*`score2oppct'^1) if riskreg == 1 & sex == 2
	replace riskage3 = min(., 51.5944 + 9.7937*ln(`score2oppct') + 0.0736*`score2oppct'^1) if riskreg == 2 & sex == 2
	replace riskage3 = min(., 46.7572 + 10.0149*ln(`score2oppct') + 0.0010*`score2oppct'^2) if riskreg == 3 & sex == 2
	replace riskage3 = min(., 18.1417 + 15.9391*ln(`score2oppct') + 0.0015*`score2oppct'^2) if riskreg == 4 & sex == 2
	replace riskage4 = min(., 41.1986 + 13.0918*ln(`score2oppct') + 0.0889*`score2oppct'^1) if riskreg == 1 & sex == 1
	replace riskage4 = min(., 39.2632 + 12.4757*ln(`score2oppct') + 0.0942*`score2oppct'^1) if riskreg == 2 & sex == 1
	replace riskage4 = min(., 33.6547 + 13.7977*ln(`score2oppct') + 0.0014*`score2oppct'^2) if riskreg == 3 & sex == 1
	replace riskage4 = min(., -6.7154 + 22.9311*ln(`score2oppct') + 0.0020*`score2oppct'^2) if riskreg == 4 & sex == 1
	replace riskage4 = min(., 49.3662 + 11.3030*ln(`score2oppct') + 0.0755*`score2oppct'^1) if riskreg == 1 & sex == 2
	replace riskage4 = min(., 48.8975 + 10.3311*ln(`score2oppct') + 0.0782*`score2oppct'^1) if riskreg == 2 & sex == 2
	replace riskage4 = min(., 43.7821 + 10.5802*ln(`score2oppct') + 0.0011*`score2oppct'^2) if riskreg == 3 & sex == 2
	replace riskage4 = min(., 13.5288 + 16.8521*ln(`score2oppct') + 0.0016*`score2oppct'^2) if riskreg == 4 & sex == 2
	replace riskage5 = min(., 43.0595 + 12.6198*ln(`score2oppct') + 0.0859*`score2oppct'^1) if riskreg == 1 & sex == 1
	replace riskage5 = min(., 41.1908 + 12.0276*ln(`score2oppct') + 0.0910*`score2oppct'^1) if riskreg == 2 & sex == 1
	replace riskage5 = min(., 35.8089 + 13.2877*ln(`score2oppct') + 0.0013*`score2oppct'^2) if riskreg == 3 & sex == 1
	replace riskage5 = min(., -3.0478 + 22.0718*ln(`score2oppct') + 0.0020*`score2oppct'^2) if riskreg == 4 & sex == 1
	replace riskage5 = min(., 50.5126 + 11.0270*ln(`score2oppct') + 0.0736*`score2oppct'^1) if riskreg == 1 & sex == 2
	replace riskage5 = min(., 50.0544 + 10.0797*ln(`score2oppct') + 0.0763*`score2oppct'^1) if riskreg == 2 & sex == 2
	replace riskage5 = min(., 45.0705 + 10.3163*ln(`score2oppct') + 0.0010*`score2oppct'^2) if riskreg == 3 & sex == 2
	replace riskage5 = min(., 15.5817 + 16.4259*ln(`score2oppct') + 0.0015*`score2oppct'^2) if riskreg == 4 & sex == 2
	replace riskage6 = min(., 44.7905 + 12.1806*ln(`score2oppct') + 0.0830*`score2oppct'^1) if riskreg == 1 & sex == 1
	replace riskage6 = min(., 42.9841 + 11.6105*ln(`score2oppct') + 0.0881*`score2oppct'^1) if riskreg == 2 & sex == 1
	replace riskage6 = min(., 37.8105 + 12.8141*ln(`score2oppct') + 0.0013*`score2oppct'^2) if riskreg == 3 & sex == 1
	replace riskage6 = min(., 0.3585 + 21.2742*ln(`score2oppct') + 0.0019*`score2oppct'^2) if riskreg == 4 & sex == 1
	replace riskage6 = min(., 51.6042 + 10.7642*ln(`score2oppct') + 0.0719*`score2oppct'^1) if riskreg == 1 & sex == 2
	replace riskage6 = min(., 51.1562 + 9.8403*ln(`score2oppct') + 0.0745*`score2oppct'^1) if riskreg == 2 & sex == 2
	replace riskage6 = min(., 46.2966 + 10.0652*ln(`score2oppct') + 0.0010*`score2oppct'^2) if riskreg == 3 & sex == 2
	replace riskage6 = min(., 17.5345 + 16.0207*ln(`score2oppct') + 0.0015*`score2oppct'^2) if riskreg == 4 & sex == 2
	foreach i of numlist 1(1)6 {
		replace riskage`i' = 0 if riskage`i' < 0
	}
end


exit		/* force exit to ignore example code below */

* examples
* CVD incidence risk calculations
score2oprisk
score2oprisk, riskreg(1) replace
score2oprisk, percent
score2oprisk, riskreg(1) replace percent

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
score2oprisk
score2oprisk, riskreg(1) replace
score2oprisk, riskreg(2) replace

* managing generated variables
des _SCORE2_START - _SCORE2_END
drop _SCORE2_START - _SCORE2_END
