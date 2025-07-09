*! Version : 0.00 12Sep2021
*! Authors : Dr Stephen Kaptoge
*!				 Dr Lisa Pennells
*! Address : Cardiovascular Epidemiology Unit
*!           Dept of Public Health and Primary Care
*!           University of Cambridge
*!           Strangeways Research Laboratory
*!           Worts Causeway
*!           Cambridge CB1 8RN
*!           UK

/*----------------------------------------------------------------------------------------------------------
S Kaptoge: Jul 2020
* Stata program for calculation of SCORE1 CVD risk score, with region specific recalibration.

* required variables: 
tchol = total cholesterol in mmol/l
sbp = systolic blood pressure in mmHg
ages = baseline age in years
sex = sex (1 = male, 2 = female)
smallbin = smoking as a binary variable (current = 1, other = 0)
hxdiabbin = history of diabetes as a binary variable (definite = 1, other = 0)

score1rreg =  the risk region (0 = low, 1 = high) 
that applies to the country of residence of the individuals in the dataset 
(see separate table for risk region classification of each country)

relevant descriptive notes.
Jul 2015:	Inital script coded as stata dofile for risk score calculation
Sep 2021:	Above converted to current wrapper program for general use
Oct 2021:	Added optional calculation of risk age
Mar 2022:	Added option multdiab(integer 1) to optionally include/exclude consideration of diabetes status
Apr 2022:	Added [if] option to be applied at the end of calculations
------------------------------------------------------------------------------------------------------------*/

capture program drop score1risk
program define score1risk, rclass
	
	version 13
	syntax [if] [, CCODEs score1rreg(numlist integer max=1 >=1 <=2) replace percent shortlab rounded(integer 0) ///
	dppscale(real 0) noWEB WEBVERSion(integer 2) DETAILs RISKAGE MULTDIAB(integer 1)]
	
	marksample touse, novarlist
	
	* preserve data in case calculation fails
	preserve
	* web settings
	local webopt = cond("`web'" == "noweb", 0, `webversion')
	local quinoi = "qui"
	`quinoi' {		/* start of quinoi */
		* drop variables from previous calculations if they exist
		capture ds _SCORE1_START - _SCORE1_END
		if _rc == 0 {
			drop `r(varlist)'
		}
		
		* check if required predictor variables exist
		local expvars = "sex ages sbp tchol hdl smallbin hxdiabbin"
		noi confirm numeric variable `expvars'
		capture noi confirm variable score1rreg, exact
		if _rc~=0 {
			if missing("`score1rreg'") {
				noi di _newline as text "Please enter SCORE1 risk region (1 = Low, 2 = High)"
				qui di _request(_answer)
				local score1rreg = `answer'
			}
			gen score1rreg = `score1rreg'
			label define score1rreg 1 "Low" 2 "High", modify
			label define score1rreg 2 "Low risk region" 2 "High risk region", modify
			label values score1rreg score1rreg
			label variable score1rreg "SCORE1 risk region"
			local orderreg = "order score1rreg, after(_SCORE1_START)"
		}
		else {
			capture confirm numeric variable score1rreg, exact
			if _rc == 0 {
				 if "`replace'" == "replace" & !missing("`score1rreg'") { 
					replace score1rreg = `score1rreg'
				}
			}
			else {
				di as error "variable score1rreg is not numeric"
				des score1rreg
				exit 197
			}
		}
				
		* generate marker for start of SCORE1 added variables
		gen _SCORE1_START = .
		`orderreg'		/* order risk region variable if generated above */
		
		* drop existing derived variables
		local dropvars1 = ""
		local dropvars2 = "score03_tchol_*"
		local dropvars3 = "score1"
		local dropvars = "`dropvars1' `dropvars2' `dropvars3'"
		local varabbrev = c(varabbrev)
		nobreak {
			set varabbrev off
			foreach var in `dropvars' {
				capture drop `var'
			}
			set varabbrev `varabbrev'
		}
		
		* Define risk level for countries either High risk or Low risk.
		capture confirm string variable country
		if _rc == 0 {
			capture drop risk_level
			gen risk_level = .
			local countries01 = "Andorra Austria Belgium Cyprus Denmark Finland France Germany Greece Iceland Ireland Israel Italy"
			local countries02 = "Luxembourg Malta Monaco Netherlands Norway Portugal San*Marino Slovenia Spain Sweden Switzerland"
			local countries03 = "United*Kingdom Northern*Ireland Scotland UK"
			local countries04 = "France Italy Spain Japan Greek Greece Belgium Netherland"
			local countries0 = "`countries01' `countries02' `countries03' `countries04'"
			foreach cou in `countries0' {
				replace risk_level = 0 if strmatch(upper(country), upper("*`cou'*"))
			}
			replace risk_level = 1 if missing(risk_level)
			table country, by(risk_level) concise
		}
		
		capture drop score03_tchol_wchd
		gen score03_tchol_wchd = .
		replace score03_tchol_wchd = 0.71*smallbin + 0.24*(tchol - 6) + 0.018*(sbp - 120)                /* weight for CHD */
		capture drop score03_tchol_s0age0_chd 
		gen score03_tchol_s0age0_chd = .
		replace score03_tchol_s0age0_chd = exp(-exp(-22.1)*((ages-20)^4.71)) if sex == 1 & score1rreg == 1 
		replace score03_tchol_s0age0_chd = exp(-exp(-29.8)*((ages-20)^6.36)) if sex == 2 & score1rreg == 1
		replace score03_tchol_s0age0_chd = exp(-exp(-21.0)*((ages-20)^4.62)) if sex == 1 & score1rreg == 2 
		replace score03_tchol_s0age0_chd = exp(-exp(-28.7)*((ages-20)^6.23)) if sex == 2 & score1rreg == 2 
		capture drop score03_tchol_s0age1_chd 
		gen score03_tchol_s0age1_chd = .
		replace score03_tchol_s0age1_chd = exp(-exp(-22.1)*((ages-10)^4.71)) if sex == 1 & score1rreg == 1 
		replace score03_tchol_s0age1_chd = exp(-exp(-29.8)*((ages-10)^6.36)) if sex == 2 & score1rreg == 1
		replace score03_tchol_s0age1_chd = exp(-exp(-21.0)*((ages-10)^4.62)) if sex == 1 & score1rreg == 2 
		replace score03_tchol_s0age1_chd = exp(-exp(-28.7)*((ages-10)^6.23)) if sex == 2 & score1rreg == 2 

		capture drop score03_tchol_wnonchd
		gen score03_tchol_wnonchd = .
		replace score03_tchol_wnonchd = 0.63*smallbin + 0.02*(tchol - 6) + 0.022*(sbp - 120)           /* weight for Non-CHD */
		capture drop score03_tchol_s0age0_nonchd
		gen score03_tchol_s0age0_nonchd = .
		replace score03_tchol_s0age0_nonchd = exp(-exp(-26.7)*((ages-20)^5.64)) if sex == 1 & score1rreg == 1 
		replace score03_tchol_s0age0_nonchd = exp(-exp(-31.0)*((ages-20)^6.62)) if sex == 2 & score1rreg == 1
		replace score03_tchol_s0age0_nonchd = exp(-exp(-25.7)*((ages-20)^5.47)) if sex == 1 & score1rreg == 2 
		replace score03_tchol_s0age0_nonchd = exp(-exp(-30.0)*((ages-20)^6.42)) if sex == 2 & score1rreg == 2 
		capture drop score03_tchol_s0age1_nonchd
		gen score03_tchol_s0age1_nonchd = .
		replace score03_tchol_s0age1_nonchd = exp(-exp(-26.7)*((ages-10)^5.64)) if sex == 1 & score1rreg == 1 
		replace score03_tchol_s0age1_nonchd = exp(-exp(-31.0)*((ages-10)^6.62)) if sex == 2 & score1rreg == 1
		replace score03_tchol_s0age1_nonchd = exp(-exp(-25.7)*((ages-10)^5.47)) if sex == 1 & score1rreg == 2 
		replace score03_tchol_s0age1_nonchd = exp(-exp(-30.0)*((ages-10)^6.42)) if sex == 2 & score1rreg == 2 

		capture drop score03_tchol_sage?_* 
		gen score03_tchol_sage0_chd = score03_tchol_s0age0_chd^exp(score03_tchol_wchd)
		gen score03_tchol_sage1_chd = score03_tchol_s0age1_chd^exp(score03_tchol_wchd)

		gen score03_tchol_sage0_nonchd = score03_tchol_s0age0_nonchd^exp(score03_tchol_wnonchd)
		gen score03_tchol_sage1_nonchd = score03_tchol_s0age1_nonchd^exp(score03_tchol_wnonchd)

		capture drop score03_tchol_?10age_*
		gen score03_tchol_s10age_chd = score03_tchol_sage1_chd/score03_tchol_sage0_chd
		gen score03_tchol_s10age_nonchd = score03_tchol_sage1_nonchd/score03_tchol_sage0_nonchd

		gen score03_tchol_p10age_chd = 1 - score03_tchol_s10age_chd
		gen score03_tchol_p10age_nonchd = 1 - score03_tchol_s10age_nonchd

		capture drop score03_tchol
		gen score03_tchol = .
		*replace score03_tchol = (score03_tchol_p10age_chd + score03_tchol_p10age_nonchd)
		replace score03_tchol = (1 - score03_tchol_s10age_chd*score03_tchol_s10age_nonchd)

		capture drop s010_score03_tchol*
		gen s010_score03_tchol_chd = exp(ln(score03_tchol_s10age_chd)/exp(score03_tchol_wchd))
		gen s010_score03_tchol_nonchd = exp(ln(score03_tchol_s10age_nonchd)/exp(score03_tchol_wnonchd))
		gen s010_score03_tchol = s010_score03_tchol_chd*s010_score03_tchol_nonchd 
		capture drop cent_score03_tchol_w*chd
		gen cent_score03_tchol_wchd = 0
		gen cent_score03_tchol_wnonchd = 0
		order s010_score03_tchol* cent_score03_tchol_w*chd, before(score03_tchol)

		*If you consider the status of history of diabetes, the paper suggest to calculate the final score by two times for male with diabetes and four times
		* for female with diabetes.
		if `multdiab' > 0 {
			*replace score03_tchol = score03_tchol*2 if sex == 1 & hxdiabbin == 1
			*replace score03_tchol = score03_tchol*4 if sex == 2 & hxdiabbin == 1
		}
		* SK comment 2013: the above statement applies to the hazard rates so 1 - survival^HR should give better approximation than the direct multiplication)
		if `multdiab' > 0 {
			replace score03_tchol = (1 - (1-(score03_tchol))^2) if sex == 1 & hxdiabbin == 1
			replace score03_tchol = (1 - (1-(score03_tchol))^4) if sex == 2 & hxdiabbin == 1
		}
		label variable score03_tchol "SCORE 2003 10-yr CVD mortality risk"

		* clone final variables for SCORE1 risk
		clonevar score1 = score03_tchol
		
		* label variables informatively
		label variable score1 "SCORE1 10-yr CVD mortality risk"
		if "`shortlab'" == "shortlab" {
			label variable score1 "SCORE1-CVDm"
		}
		
		if "`percent'" == "percent" {
			local riskvars = "score1"
			foreach var of varlist `riskvars' {
				local varlab: variable label `var'
				replace `var' = 100*`var'		/* convert to percent risk */
				label variable `var' "`varlab' (%)"
			}
		}
		
		* calculate risk age if requested
		if !missing("`riskage'") {
			score1riskage if !missing(score1)		/* call wrapper program below to calculate risk age */
		}

		* generate marker for end of SCORE1 added variables
		gen _SCORE1_END = .
	}		/* end of quinoi */
	
	* cancel restore if successful
	restore, not
	
	* keep minimal variables in dataset if nodetails
	if inlist(`webopt', 1, 2) {
		qui ds _SCORE1_START - _SCORE1_END
		local dropvars = r(varlist)
		local calvars = cond(inlist(`webopt', 1), "score1", "score1")
		if !missing("`riskage'") {
			local calvars = "`calvars' ages riskage?"
		}
		if "`details'" ~= "details" {
			local keepvars = "_SCORE1_START score1rreg `calvars' _SCORE1_END"
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
	qui recode _SCORE1_START - _SCORE1_END (* = .) if !(`touse')
	des _SCORE1_START - _SCORE1_END
	preserve
	qui keep if !missing(score1)
	local calvars = cond(inlist(`webopt', 0, 1), "score1", "score1")
	if !missing("`riskage'") {
		local calvars = "`calvars' ages riskage?"
	}
	bysort score1rreg: tabstat `calvars', by(sex) stats(n mean sd min p25 p50 p75 max) ///
	col(stats) longstub varw(16)
	if !missing("`riskage'") {
		notes `calvars'
	}
	restore
	
end

* wrapper program to calculate risk age for SCORE1
capture program drop score1riskage
program define score1riskage, rclass
	
	version 13
	syntax [if] [in] [, *]
	marksample touse, novarlist
	
	local dropvars = "riskage1 riskage2 riskage3 riskage4 riskage5 riskage6"
	foreach var in `dropvars' {
		capture drop `var'
	}
	foreach i of numlist 1(1)6 {
		gen riskage`i' = .
		label variable riskage`i' "SCORE1 risk age`i'"
	}
	notes riskage1: "SCORE1 risk age1: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage2: "SCORE1 risk age2: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage3: "SCORE1 risk age3: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.4 (TCHOL = 4.7[M] 4.9[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage4: "SCORE1 risk age4: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 135, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage5: "SCORE1 risk age5: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 130, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	notes riskage6: "SCORE1 risk age6: Ideal RFx (smallbin = 0, hxdiabbin = 0, SBP = 125, nonHDL = 3.9 (TCHOL = 5.2[M] 5.4[F] & HDL = 1.3[M] 1.5[F]))"
	tempvar score1pct
	summ score1
	gen `score1pct' = cond(r(max) < 1, score1*100, score1)
	* risk age estimation equations derived separately (see score1risk_test.do)
	replace riskage1 = min(., 42.2822 + 3.9597*ln(`score1pct') + 8.5205*`score1pct'^.5) if score1rreg == 1 & sex == 1
	replace riskage1 = min(., 37.8333 + 3.6680*ln(`score1pct') + 7.5700*`score1pct'^.5) if score1rreg == 2 & sex == 1
	replace riskage1 = min(., 51.1122 + 3.6477*ln(`score1pct') + 6.8470*`score1pct'^.5) if score1rreg == 1 & sex == 2
	replace riskage1 = min(., 48.7457 + 3.6244*ln(`score1pct') + 6.3015*`score1pct'^.5) if score1rreg == 2 & sex == 2
	replace riskage2 = min(., 42.7570 + 3.9842*ln(`score1pct') + 8.8658*`score1pct'^.5) if score1rreg == 1 & sex == 1
	replace riskage2 = min(., 38.3354 + 3.7119*ln(`score1pct') + 7.8209*`score1pct'^.5) if score1rreg == 2 & sex == 1
	replace riskage2 = min(., 51.5656 + 3.6630*ln(`score1pct') + 7.1240*`score1pct'^.5) if score1rreg == 1 & sex == 2
	replace riskage2 = min(., 49.2213 + 3.6473*ln(`score1pct') + 6.5263*`score1pct'^.5) if score1rreg == 2 & sex == 2
	replace riskage3 = min(., 43.2262 + 4.0070*ln(`score1pct') + 9.2332*`score1pct'^.5) if score1rreg == 1 & sex == 1
	replace riskage3 = min(., 38.8289 + 3.7515*ln(`score1pct') + 8.0912*`score1pct'^.5) if score1rreg == 2 & sex == 1
	replace riskage3 = min(., 61.7188 + 8.8597*ln(`score1pct') + 0.5163*ln(`score1pct')) if score1rreg == 1 & sex == 2
	replace riskage3 = min(., 49.6777 + 3.6632*ln(`score1pct') + 6.7738*`score1pct'^.5) if score1rreg == 2 & sex == 2
	replace riskage4 = min(., 41.7227 + 3.9291*ln(`score1pct') + 8.3587*`score1pct'^.5) if score1rreg == 1 & sex == 1
	replace riskage4 = min(., 37.2099 + 3.6103*ln(`score1pct') + 7.4260*`score1pct'^.5) if score1rreg == 2 & sex == 1
	replace riskage4 = min(., 50.8132 + 3.6425*ln(`score1pct') + 6.6977*`score1pct'^.5) if score1rreg == 1 & sex == 2
	replace riskage4 = min(., 48.3747 + 3.6085*ln(`score1pct') + 6.1560*`score1pct'^.5) if score1rreg == 2 & sex == 2
	replace riskage5 = min(., 42.1987 + 3.9551*ln(`score1pct') + 8.6897*`score1pct'^.5) if score1rreg == 1 & sex == 1
	replace riskage5 = min(., 37.7168 + 3.6581*ln(`score1pct') + 7.6612*`score1pct'^.5) if score1rreg == 2 & sex == 1
	replace riskage5 = min(., 51.2431 + 3.6493*ln(`score1pct') + 6.9761*`score1pct'^.5) if score1rreg == 1 & sex == 2
	replace riskage5 = min(., 48.8476 + 3.6304*ln(`score1pct') + 6.3714*`score1pct'^.5) if score1rreg == 2 & sex == 2
	replace riskage6 = min(., 42.6688 + 3.9791*ln(`score1pct') + 9.0422*`score1pct'^.5) if score1rreg == 1 & sex == 1
	replace riskage6 = min(., 38.2147 + 3.7013*ln(`score1pct') + 7.9154*`score1pct'^.5) if score1rreg == 2 & sex == 1
	replace riskage6 = min(., 51.7112 + 3.6709*ln(`score1pct') + 7.2487*`score1pct'^.5) if score1rreg == 1 & sex == 2
	replace riskage6 = min(., 49.3113 + 3.6497*ln(`score1pct') + 6.6034*`score1pct'^.5) if score1rreg == 2 & sex == 2
	foreach i of numlist 1(1)6 {
		replace riskage`i' = 0 if riskage`i' < 0
	}
end


exit		/* force exit to ignore example code below */

* examples
* CVD incidence risk calculations
score1risk
score1risk, score1rreg(1) replace
score1risk, percent
score1risk, score1rreg(1) replace percent

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
score1risk
score1risk, score1rreg(1) replace
score1risk, score1rreg(2) replace

* managing generated variables
des _SCORE1_START - _SCORE1_END
drop _SCORE1_START - _SCORE1_END
