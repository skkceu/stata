*! Version : 0.00 12Dec2021
*! Author  : Dr Stephen Kaptoge
*! Address : Cardiovascular Epidemiology Unit
*!           Dept of Public Health and Primary Care
*!           University of Cambridge
*!           Strangeways Research Laboratory
*!           Worts Causeway
*!           Cambridge CB1 8RN
*!           UK

/*---------------------------------------------------------------------------------------------------------------
S Kaptoge Dec 2021: wrapper program to apply eGFR/CKDACR patches for SCORE2 or SCORE2-OP optionally
-----------------------------------------------------------------------------------------------------------------*/

* wrapper program to apply eGFR/CKDACR patches for SCORE2 or SCORE2-OP optionally
capture program drop score2egfr
program define score2egfr, rclass
	version 11
	syntax [, varpfix(name) riskreg(numlist integer max=1 >=1 <=104) replace percent ///
	egfrvar(varname numeric) acrvar(varname numeric) noscore2risk details]
	
	marksample touse, novarlist
	
	* define defaults and drop variables
	local varpfix = cond(missing("`varpfix'"), "score2", "`varpfix'")
	local varpfix = cond(!strmatch("`varpfix'", "*_"), "`varpfix'_", "`varpfix'")
	local egfrvar = cond(missing("`egfrvar'"), "egfr_ckdepi", "`egfrvar'")
	local acrvar = cond(missing("`acrvar'"), "acr", "`acrvar'")
	local marker_start = upper("_`varpfix'") + upper("_egfr_patch_start")
	local marker_end = upper("_`varpfix'") + upper("_egfr_patch_end")
	* preserve data in case calculation fails
	preserve
	* calculate score 2
	if "`score2risk'" ~= "noscore2risk" {
		score2risk, riskreg(`riskreg') detail score2op `replace'
	}
	capture drop `marker_start' - `marker_end'
	gen `marker_start' = .
	local dropvars1 = "scale1 scale2 egfr_pred acr_pred log8acr log8acr_pred"
	local dropvars2 = "`varpfix'ucal `varpfix'egfr_* `varpfix'ckdacr_* `varpfix'ucal_*"
	local dropvars3 = "cal?_`varpfix'egfr_* cal?_`varpfix'ckdacr_* cal?_`varpfix'ucal_*"
	local dropvars4 = "pct_`varpfix'egfr_* pct_`varpfix'ckdacr_* pct_`varpfix'ucal_*"
	local dropvars5 = "pct_cal?_`varpfix'egfr_* pct_cal?_`varpfix'ckdacr_* pct_cal?_`varpfix'ucal_*"
	local dropvars = "`dropvars1' `dropvars2' `dropvars3' `dropvars4' `dropvars5'"
	foreach var in `dropvars' {
		capture drop `var'
	}
	local quinoi = "noi"
	`quinoi' {
		local ucalvar = "`varpfix'ucal"
		local ucallab = cond(inlist("`varpfix'", "score2_"), "ERFC", "CONOR")
		clonevar `ucalvar' = `varpfix'cal0_m1_risk
		local scorevar = subinstr("`varpfix'", "_", "", 1)
		local scorelab = upper("`scorevar'")
		local scale1 = "`varpfix'scale1"
		local scale2 = "`varpfix'scale2"
		capture drop scale1 scale2
		clonevar scale1 = `scale1'
		clonevar scale2 = `scale2'
		* ERFC model recalibrated SCORE2
		local predrisk = "`ucalvar'"
		label variable `predrisk' "`ucallab'"
	*	recalib, prerisk(`predrisk') genvar(`scorevar') labelpre("`ucallab'") labelcal("`scorelab'")
		* expected eGFR	
		gen egfr_pred = 87.8980 - 3.7891*(ages - 60)/5 - 0.7023*(sex == 2) - 0.2941*(tchol - 6) + 1.0960*(hdl - 1.3)/0.5 - 0.1364*(sbp - 120)/20 + 0.1205*(hxdiabbin == 1) + 1.3211*(smallbin == 1) + 0.0555*(ages - 60)/5*(tchol - 6) + 0.1717*(ages - 60)/5*(hdl - 1.3)/0.5 + 0.0059*(ages - 60)/5*(sbp - 120)/20 - 0.8994*(ages - 60)/5*(hxdiabbin == 1) + 0.2181*(ages - 60)/5*(smallbin == 1)
		label variable egfr_pred "Predicted eGFR (ml/min/1.73m^2)"
		local egfrobs = "`egfrvar'"
		local egfrexp = "egfr_pred"
		* eGFR patch to recalibrated SCORE2
		local predrisk = "`scorevar'"
		local patchvar = "`predrisk'_egfr_p"
		local patcheqn_score2   = "1 - (1 - \`predrisk')^exp(0.4713*(min(`egfrobs', 60)/-15 - min(`egfrexp', 60)/-15) + 0.0956*(min(max(`egfrobs' - 60, 0), 30)/-15 - min(max(`egfrexp' - 60, 0), 30)/-15) - 0.3308*(max(`egfobs' - 90, 0)/-15 - max(`egfrexp' - 90, 0)/-15) - 0.0802*(ages - 60)/5*(min(`egfrobs', 60)/-15 - min(`egfrexp', 60)/-15) + 0.0088*(ages - 60)/5*(min(max(`egfobs' - 60, 0), 30)/-15 - min(max(`egfrexp' - 60, 0), 30)/-15) - 0.0224*(ages - 60)/5*(max(`egfrobs' - 90, 0)/-15 - max(`egfrexp' - 90, 0)/-15)) if !missing(`egfrobs')"
		local patcheqn_score2op = "1 - (1 - \`predrisk')^exp(0.3072*(min(`egfrobs', 60)/-15 - min(`egfrexp', 60)/-15) + 0.0942*(min(max(`egfrobs' - 60, 0), 30)/-15 - min(max(`egfrexp' - 60, 0), 30)/-15) - 0.4616*(max(`egfrobs' - 90, 0)/-15 - max(`egfrexp' - 90, 0)/-15) - 0.0127*(ages - 73)*(min(`egfrobs', 60)/-15 - min(`egfrexp', 60)/-15) - 0.0098*(ages - 73)*(min(max(`egfrobs' - 60, 0), 30)/-15 - min(max(`egfrexp' - 60, 0), 30)/-15) - 0.0075*(ages - 73)*(max(`egfrobs' - 90, 0)/-15 - max(`egfrexp' - 90, 0)/-15)) if !missing(`egfrobs')"
		local patcheqn = "\`patcheqn_`scorevar''"
		gen `patchvar' = `patcheqn'
		label variable `patchvar' "`scorelab'+eGFR_p"
		recalib, prerisk(`patchvar') genvar(`scorevar'_egfr_v01) norescale labelpre("`scorelab'+eGFR_p") labelcal("`scorelab'+eGFR_p-RC")
		* eGFR patch to uncalibrated SCORE2
		local predrisk = "`ucalvar'"
		local patchvar = "`predrisk'_egfr_p"
		gen `patchvar' = `patcheqn'
		label variable `patchvar' "`ucallab'+eGFR_p"
		recalib, prerisk(`patchvar') genvar(`scorevar'_egfr_v02) labelpre("`ucallab'+eGFR_p") labelcal("`ucallab'+eGFR_p+RC")
		* expected ACR
		gen acr_pred = 8^(1 - 0.0225 + 0.0159*(ages - 60)/5 + 0.0308*(sex == 2) + 0.0185*(tchol - 6) - 0.0274*(hdl - 1.3)/0.5 + 0.1339*(sbp - 120)/20 + 0.2171*(hxdiabbin == 1) + 0.0629*(smallbin == 1) - 0.0062*(ages - 60)/5*(tchol - 6) + 0.0003*(ages - 60)/5*(hdl - 1.3)/0.5 + 0.0008*(ages - 60)/5*(sbp - 120)/20 - 0. 0109*(ages - 60)/5*(hxdiabbin == 1) + 0.0085*(ages - 60)/5*(smallbin == 1) + 0.4057*min(`egfrobs' - 60, 0)/-15 + 0.0597*min(max(`egfrobs' - 60, 0), 30)/-15 - 0.0916*max(`egfrobs' - 90, 0)/-15) if !missing(`egfrobs')
		label variable acr_pred "Predicted ACR (mg/g)"
		gen log8acr_pred = ln(acr_pred)/ln(8)
		label variable log8acr_pred "Log8 Predicted ACR (mg/g)"
		* CKDACR patch to eGFR patch
		local acrobs = "`acrvar'"
		local acrexp = "acr_pred"
		capture noi confirm variable `acrobs', exact
		if _rc == 0 {
			local log8acrvar = "log8`acrobs'"
			gen `log8acrvar' = ln(`acrobs')/ln(8)
			label variable `log8acrvar' "Log8 `: variable label `acrobs''"
			* CKDACR patch to calibrated SCORE2 eGFR patch
			local acrobs = "8.85*`acrobs'"		/* convert to mg/g as used in acr_pred formulae */
			local predrisk = "`scorevar'_egfr_v01"
			local patchvar = "`scorevar'_ckdacr_p_v01"
			local patcheqn_score2   = "1 - (1 - \`predrisk')^exp(0.2432*(ln(`acrobs')/ln(8) - ln(`acrexp')/ln(8))) if !missing(`acrobs')"
			local patcheqn_score2op = "1 - (1 - \`predrisk')^exp(0.2370*(ln(`acrobs')/ln(8) - ln(`acrexp')/ln(8))) if !missing(`acrobs')"
			local patcheqn = "\`patcheqn_`scorevar''"
			gen `patchvar' = `patcheqn'
			label variable `patchvar' "`scorelab'+eGFR_p-RC+CKDACR_p"
			recalib, prerisk(`patchvar') genvar(`scorevar'_ckdacr_v01) norescale ///
			labelpre("`scorelab'+eGFR_p-RC+CKDACR_p") labelcal("`scorelab'+eGFR_p-RC+CKDACR_p-RC")
			local predrisk = "`scorevar'_egfr_v02"
			local patchvar = "`scorevar'_ckdacr_p_v02"
			gen `patchvar' = `patcheqn'
			label variable `patchvar' "`ucallab'+eGFR_p+RC+CKDACR_p"
			recalib, prerisk(`patchvar') genvar(`scorevar'_ckdacr_v02) norescale ///
			labelpre("`ucallab'+eGFR_p+RC+CKDACR_p") labelcal("`ucallab'+eGFR_p+RC+CKDACR_p-RC")
			* CKDACR patch to uncalibrated SCORE2 eGFR patch
			local predrisk = "`scorevar'_ucal_egfr_p"
			local patchvar = "`scorevar'_ucal_ckdacr_p"
			gen `patchvar' = `patcheqn'
			label variable `patchvar' "`ucallab'+eGFR_p+RC+CKDACR_p"
			recalib, prerisk(`patchvar') genvar(`scorevar'_ckdacr_v03) ///
			labelpre("`ucallab'+eGFR_p+CKDACR_p") labelcal("`ucallab'+eGFR_p+CKDACR_p+RC")
		}
	}
	gen `marker_end' = .
	des `marker_start' - `marker_end'
	* cancel restore if successful
	restore, not
	
	if "`percent'" == "percent" {
		local riskvars = "`varpfix'ucal score2 score2?? `varpfix'*_v0?"
		foreach var of varlist `riskvars' {
			local varlab: variable label `var'
			replace `var' = 100*`var'		/* convert to percent risk */
			label variable `var' "`varlab' (%)"
		}
	}
	
	* describe and summarize risk variables
	qui recode `marker_start' - `marker_end' (* = .) if !(`touse')
	des `marker_start' - `marker_end'
	preserve
	qui keep if !missing(score2)
	local calvars = cond(inlist(`webopt', 0, 1), "`varpfix'ucal score2 score2?? `varpfix'*_v0?", "`varpfix'ucal score2 score2?? `varpfix'*_v0?")
	if !missing("`riskage'") {
		local calvars = "`calvars' ages riskage?"
	}
	des `calvars'
	capture ds riskreg, exact
	if _rc ~= 0 {
		capture ds score2_riskreg
	}
	local byriskreg = cond(_rc == 0, "bysort `r(varlist)': ", "")
	`byriskreg' tabstat `calvars', by(sex) stats(n mean sd min p25 p50 p75 max) ///
	col(stats) longstub varw(16)
	if !missing("`riskage'") {
		notes `calvars'
	}
	restore

end

* small program to apply recalibration parameters to predicted risks
capture program drop recalib
program define recalib
	version 11
	syntax [varlist(default=none)] [if] [in] [, prerisk(varname numeric) ///
	genvar(namelist max=1) labelpre(string) labelcal(string) addlabel(string) ///
	norescale replace b_cons(varname) b_slope(varname)]
	tempvar loglogpre
	local prerisk = cond(missing("`prerisk'"), "p10_m1_ep_cv", "`prerisk'")
	local genvar = cond(missing("`genvar'"), "score2", "`genvar'")
	if missing("`b_cons'") {
		local b_cons = cond("`rescale'" == "norescale", "0", "scale1")
	}
	if missing("`b_slope'") {
		local b_slope = cond("`rescale'" == "norescale", "1", "scale2")
	}
	foreach i of numlist 2(1)2 {
		local cpfix = "cal`i'_"		/* 1 = country-specific, 2 = ESC regions */
		local calrisk = "`cpfix'`prerisk'"
		gen `loglogpre' = ln(-ln(1 - `prerisk'))
		if "`replace'" == "replace" {
			capture drop pct_`prerisk'
			capture drop `calrisk'
			capture drop pct_`calrisk'
			capture drop `genvar'
			capture drop pct_`genvar'
		}
		gen `calrisk' = 1 - exp(-exp(`b_cons' + `b_slope'*`loglogpre'))
		drop `loglogpre'
		local labelpre = cond(missing("`labelpre'"), "ERFC", "`labelpre'")
		local labelcal = cond(missing("`labelcal'"), "SCORE2", "`labelcal'")
		label variable `prerisk' "`labelpre'`addlabel'"
		label variable `calrisk' "`labelcal'`addlabel'"
		clonevar pct_`prerisk' = `prerisk'
		clonevar pct_`calrisk' = `calrisk'
		replace pct_`prerisk' = 100*`prerisk'
		replace pct_`calrisk' = 100*`calrisk'
		clonevar `genvar' = `calrisk'
		clonevar pct_`genvar' = pct_`calrisk'
	}
end

exit		/* force exit to ignore code below */

* example calls
recalib
recalib, prerisk(p10_m1_ep_cv) genvar(score2)
recalib, prerisk(p10_m1_ep_cv_egfr) genvar(score2_egfr) addlabel("+eGFR")

exit		/* force exit to ignore code below */

* example call
score2egfr
score2egfr, varpfix(score2)
score2egfr, varpfix(score2op)
score2egfr, riskreg(1) replace
score2egfr, riskreg(1) replace varpfix(score2)
score2egfr, riskreg(1) replace varpfix(score2op)

score2risk, riskreg(1) replace detail score2op
score2egfr, noscore2risk
score2egfr, noscore2risk varpfix(score2)
score2egfr, noscore2risk varpfix(score2op)

score2egfr, riskreg(1) replace percent
score2egfr, riskreg(1) replace percent varpfix(score2)
score2egfr, riskreg(1) replace percent varpfix(score2op)

