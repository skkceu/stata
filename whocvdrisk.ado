*! Version : 0.00 12Feb2019
*! Authors : Dr Stephen Kaptoge
*! Address : Cardiovascular Epidemiology Unit
*!           Dept of Public Health and Primary Care
*!           University of Cambridge
*!           Strangeways Research Laboratory
*!           Worts Causeway
*!           Cambridge CB1 8RN
*!           UK

/*----------------------------------------------------------------------------------------------------------
S Kaptoge: Feb 2019
* Stata program for calculation of WHO CVD risk score, with country/region specific recalibration.
relevant descriptive notes.
Jan 2017:	Inital script coded as stata dofile for risk score calculation in ERFC data
Jun 2018:	Added recalibration procedure for single year to dofile 
Dec 2018:	Added recalibration procedure for multiple years to dofile 
Feb 2019:	Above converted to current wrapper program for general use
Sep 2020:	Modified web folder structure to better suit future update of recalibration
Feb 2022:	Made efficient for webuse by skipping unecessary inputs
------------------------------------------------------------------------------------------------------------*/
capture program drop whocvdrisk
program define whocvdrisk, rclass
	
	version 13
	syntax [, CCODEs ratesref(string) ratesyrs(numlist integer) acuteinc(integer 0) ///
	dosvar(varname numeric) gbdyear(integer 2017) calyear(numlist integer max=1) attime(integer 10) percent ///
	shortlab norenpfix mergeinp(integer -1) mergekeep(numlist integer >=1 <=5) pause SENSCAL(integer 0) noWEB noEXTRA ///
	noFIXCALYEAR WEBVERSion(integer 2) CALRUN(integer 1) CALPARMS(string) CALINPUTS(string) ///
	SAVECALPARMS(string) SAVECALINPUTS(string) OTHVARS(string) DETAILs TABPARMs noRELABEL]
	
	* assume data in memory and perform risk score calculation and recalibration
	if missing("`ccodes'") {
		capture noi confirm variable ccode, exact
		if _rc ~= 0 {
			di as error "Country code variable (ccode) is required. To lookup from a list of country codes"
			di as error "and then generate ccode variable with relevant code, type"
			di as text "whocvdrisk, ccodes"
			exit 197
		}
	}
	if !inlist(`attime', 1, 5, 10) {
		di as error "invalid attime(`attime') - recalibration only possible for 5- or 10- year CVD risk"
		exit 197
	}
	preserve
	local whodir = "V:\ERFC\Analysis\METHODS\PubHealth\SCORE_WHO"
	local weburl = "http://ceu.phpc.cam.ac.uk/software/erfc"
	local quinoi = "qui"
	`quinoi' {
		local ratesref = cond(missing("`ratesref'"), "gbd", "`ratesref'")	/* reference rates source */
		local ratesyrs = cond(missing("`ratesyrs'"), "2017(1)2016 2015(5)1990", "`ratesyrs'")	/* target year for rates */
		local acuteinc = "`acuteinc'"		/* 2 = acute incidence + deaths, 1 = acute incidence, 0 = overall incidence */
		local dosvar = cond(missing("`dosvar'"), "dos", "`dosvar'")
		local j = cond(missing("`acuteinc'"), "1", "`acuteinc'")
		local outpdir1 = "output"		/* cd to relevant output directory */
		local outpdir0 = "sens1"
		local outpdir2 = "sens2"
		* web settings
		local webopt = cond("`web'" == "noweb", 0, `webversion')
		local calrunparm = "CALRUN" + cond(`calrun' < 10, "0", "") + string(`calrun')
		local t = "`attime'"
		if inlist(`webopt', 0, 1, 2) {
			if "`web'" == "noweb" {
				local webdir = "`whodir'\WEB\whocvdrisk\\`gbdyear'\\`calrunparm'"
			}
			else {
				local webdir = "`weburl'/whocvdrisk/`gbdyear'/`calrunparm'"
			}
			local sourcedir0 = "`whodir'\dofiles"
			local sourcedir1 = "`webdir'"
			local sourcedir2 = "`webdir'"
			local k = `webopt'		/* chooose from above */
			local sourcedir = "`sourcedir`k''"
			* check if OK to use online calculations
			local timeout1 = c(timeout1)
			set timeout1 1
			tempfile dofilecode
			local dofile = "`sourcedir'\erfc_score_v03.do"
			capture qui copy "`dofile'" "`dofilecode'"
			if !inlist(_rc, 0) {
				local adodir : sysdir PLUS
				local webdir = "`adodir'"+"w"
				local dofile = "`webdir'/whocvdrisk_erfc_score_v03.ado"
				capture qui copy "`dofile'" "`dofilecode'"
			} 
			set timeout1 `timeout1'
			if missing("`calparms'") {
				local calparms = "`webdir'\whocvdrisk_t`t'_calparms_`ratesref'_sens`j'"
			}
			if missing("`calinputs'") {
				local calinputs = "`webdir'\whocvdrisk_t`t'_calinputs"
			}
			if missing("`calyear'") & missing("`fixcalyear'") {
				local calyear = `gbdyear'
			}
			* check if valid gbdyear
			local gbdyrmin = 2017		
			local gbdyrmax = 2017
			if !inrange(`gbdyear', `gbdyrmin', `gbdyrmax') {
				qui numlist "`gbdyrmin'(1)`gbdyrmax'"
				di as error "invalid gbdyear(`gbdyear') - should be one of `r(numlist)'"
				exit 197
			}
		}
		* make lookup dataset of calibration parameters
		if !missing("`savecalparms'") {
			* get parameters from my default source
			tempfile master tempres
			local firstcall = 1
			foreach ratesyr of numlist `ratesyrs' {
				*cd "`outpdir`j''\\`ratesyr'"
				local calxyrdir = cond(`attime' <= 1, "CAL01YR", cond(`attime' <= 5, "CAL05YR", ""))
				local sourcedir = "`whodir'\\`calxyrdir'\\`outpdir`j''\\`ratesyr'\\`ratesref'"
				local senscalib = cond(inlist(`senscal', 1, 2, 3), "senscal`senscal'", "calib")
				local calparms = "`sourcedir'\output\who_score_`senscalib'_ep_merged_bystrata_coefs.dta"
			*	cd "`sourcedir'"			/* can comment out */
				use "`calparms'", clear
			*	keep ccode country countryid year sex cal?_m?_cons_* cal?_m?_slope_*
				duplicates drop
				isid year ccode sex
				duplicates drop
				char _dta[pstrata]
				char _dta[strata]
				compress
				save `tempres', replace
				if `firstcall' == 1 {
					save `master', replace
					local firstcall = 0
				}
				else {
					use `master', clear
					append using `tempres'
					save, replace
				}
			}
			use `master', clear
			if !missing("`savecalparms'") {
				saveold "`savecalparms'", version(13) replace
			}
		}
		else {
			* get parameters as supplied
			use "`calparms'", clear
		}
		clonevar calyear = year				/* to use instead of year, which may exist in user's dataset */
		tempfile calparms
		save `calparms'
		if "`pause'" == "pause" {
			pause on
			pause
		}
		* make lookup dataset of calibration inputs
		if `mergeinp' >= 0 {
			if !missing("`savecalinputs'") {
				tempfile master tempres
				local firstcall = 1
				foreach ratesyr of numlist `ratesyrs' {
					*cd "`outpdir`j''\\`ratesyr'"
					local calxyrdir = cond(`attime' <= 1, "CAL01YR", cond(`attime' <= 5, "CAL05YR", ""))
					local sourcedir1 = "`whodir'\\`calxyrdir'\\`outpdir`j''\\`ratesyr'\\`ratesref'"
					local ratesrfxs1 = "`sourcedir1'\output\who_score_calib_ep_merged_bystrata.dta"
					local sourcedir2 = "`whodir'\\`calxyrdir'\data\WHO\cleaned"
					local ratesrfxs2 = "`sourcedir2'\who_cvd_rates_rfx_data_gbd_2017_v04_yr_`ratesyr'_interpolated.dta"
					local k = 2						/* chooose from above */
					local sourcedir = "`sourcedir`k''"
					local ratesrfxs = "`ratesrfxs`k''"
					*	cd "`sourcedir'"			/* can comment out */
					use "`ratesrfxs'", clear
					local othvars = "`othvars' pr??_ep_* pr?_ep_* npop"
					local othvars: list uniq othvars
					foreach othvar in `othvars' {
						capture ds `othvar'
						if _rc ~= 0 {
							local othvars: list othvars - othvar
						}
					}
					local keepvars1 = "ccode country year sex agegrp npop cdr100k_ep_* ir100k_ep_* `othvars'"
					local keepvars2 = "midage ages sbp dbp hypbin bmi hxdiabbin smallbin tchol hdl nonhdl *id expand"
					local keepvars = "`keepvars1' `keepvars2'"
					foreach var in `keepvars' {
						capture ds `var'
						if _rc ~= 0 {
							local keepvars: list keepvars - var
						}
					}
					keep `keepvars'
					local pfixvars = "midage ages sbp dbp hypbin bmi hxdiabbin smallbin tchol hdl nonhdl"
					local pfixvars: list pfixvars & keepvars
					rename (`pfixvars') inp_=
					duplicates drop
					isid year ccode sex agegrp expand
					duplicates drop
					char _dta[pstrata]
					char _dta[strata]
					compress
					save `tempres', replace
					if `firstcall' == 1 {
						save `master', replace
						local firstcall = 0
					}
					else {
						use `master', clear
						append using `tempres'
						save, replace
					}
				}
				use `master', clear
				if !missing("`savecalinputs'") {
					saveold "`savecalinputs'", version(13) replace
				}
			}
			else {
				* get inputs as supplied
				use "`calinputs'", clear
			}
			clonevar calyear = year				/* to use instead of year, which may exist in user's dataset */
			clonevar calexpand = expand		/* to use instead of expand, which may exist in user's dataset */
			label copy expand calexpand
			label values calexpand calexpand
			tempfile calinputs
			save `calinputs'
			if "`pause'" == "pause" {
				pause on
				pause
			}
		}
		* display country code mappings to GBD regions
		if "`ccodes'" == "ccodes" {
			egen tag1 = tag(country) if gbdregid <= 21 & !strmatch(ccode, "???GBDR_*")
			noi tabdisp country if tag1 == 1, by(gbdregid) concise c(ccode)
			local extra = cond(missing("`extra'"), "noextra", "`extra'")
		}
		restore
		if "`extra'" == "noextra" {
			noi di as text "Done for now ..."
			exit
		}
		* drop variables from previous calculations if they exist
		capture ds _WHO_START - _WHO_END
		if _rc == 0 {
			drop `r(varlist)'
		}
		gen _WHO_START = .
		* calculate risk score and recalibrate
		if !missing("`calyear'") {
			generat calyear = `calyear'
		}
		else {
			generat calyear = year(`dosvar')		/* year of baseline assessment */
		}
		replace calyear = 5*floor(calyear/5) if calyear <= 2015		/* assume 5-year periods for updating recalibration */
		merge m:1 calyear ccode sex using `calparms', keep(master matched) gen(_mergecal)
		if inlist(`mergeinp', 0, 1) {
			preserve
			use `calinputs', clear
			keep if inlist(calexpand , `mergeinp')
			tempfile tempcalinputs
			save `tempcalinputs', replace
			restore
			gen calexpand = `mergeinp'
			local mergekeep = cond(missing("`mergekeep'"), "1 3", "`mergekeep'")
			merge m:1 calyear ccode sex agegrp calexpand using `tempcalinputs', keep(`mergekeep') gen(_mergeinp)
		}
		capture qui do "`dofilecode'"
		local riskvars_ep_chdmi = "chdr"
		local riskvars_ep_crbv  = "strr"
		local riskvars_ep_cvdr 	= "cvdr"
		local riskvars_ep_cvdd 	= "cvdd"
		local epvars0 = "ep_chdmi ep_crbv ep_cvdr ep_cvdd"
		local epvars1 = "ep_chdmi ep_crbv ep_cvdr ep_cvdd"
		local epvars2 = "ep_chdmi ep_crbv ep_cvdr"
		local epvars = "`epvars`k''"
		local calnos0 = "1(1)2"
		local calnos1 = "1(1)2"
		local calnos2 = "2(1)2"
		local calnos = "`calnos`k''"
		tempvar loglogpre
		foreach i of numlist `calnos' {
			local cpfix = "cal`i'_"		/* 1 = country-specific, 2 = GBD regions */
			foreach epvar in `epvars' {
				foreach rvar in `riskvars_`epvar'' {
					foreach j of numlist 1(1)2 {
						local prerisk = "erfc_`rvar'_m`j'"
						local b_cons = "`cpfix'm`j'_cons_`epvar'"
						local b_slope = "`cpfix'm`j'_slope_`epvar'"
						local calrisk = "`cpfix'`prerisk'"
						capture confirm variable `b_cons', exact
						if _rc ~= 0 {
							local b_cons = .
						}
						capture confirm variable `b_slope', exact
						if _rc ~= 0 {
							local b_slope = .
						}
						local prerisklab: variable label `prerisk'
						gen `loglogpre' = ln(-ln(1-`prerisk'))
						gen `calrisk' = 1 - exp(-exp(`b_cons' + `b_slope'*`loglogpre'))
						local sfix = cond(strmatch("`rvar'", "cvd?"), substr("`rvar'", -1, .), " ")
						local prerisklab: subinstr local prerisklab "WHO-CVD" "WHO-CVD`sfix'"
						local prerisklab: subinstr local prerisklab "Stroke" "STR"
					*	label variable `prerisk' "`prerisklab'"
						local pfix = "cal`i'"
						if missing("`shortlab'") {
							local calrisklab = trim("`t'-year risk `pfix'`prerisklab'")
							local calrisklab: subinstr local calrisklab "Stroke" "STR"
						*	local calrisklab: subinstr local calrisklab " Model " "M"
							local calrisklab = trim("`calrisklab'")
						}
						else {
							local calrisklab = "`pfix'`prerisklab'"
							local calrisklab: subinstr local calrisklab "Stroke" "STR"
							local calrisklab: subinstr local calrisklab " Model " "M"
							local calrisklab = trim("`calrisklab'")
						}
						label variable `calrisk' "`calrisklab'"
						drop `loglogpre'
					}
				}	
			}
			foreach j of numlist 1(1)2 {
				capture `quinoi' confirm variable `cpfix'erfc_chdr_m`j' `cpfix'erfc_strr_m`j'
				if _rc == 0 {
					capture `quinoi' clonevar erfc_cvdx_m`j' = erfc_cvdr_m`j'
					clonevar `cpfix'erfc_cvdx_m`j' = `cpfix'erfc_cvdr_m`j'
					replace `cpfix'erfc_cvdx_m`j' = 1 - (1 - `cpfix'erfc_chdr_m`j')*(1 - `cpfix'erfc_strr_m`j')
					order erfc_cvdx_m`j', after(erfc_cvdr_m`j')
					local calrisk = "`cpfix'erfc_cvdx_m`j'"
					local calrisklab: variable label `calrisk'
					local calrisklab: subinstr local calrisklab "WHO-CVDr" "WHO-CVDx"
					label variable `calrisk' "`calrisklab'"
				}
			}
		}
		if "`percent'" == "percent" {
			local riskvars = "erfc_*_m? cal?_erfc_*"
			foreach var of varlist `riskvars' {
				local varlab: variable label `var'
				replace `var' = 100*`var'		/* convert to percent risk */
				label variable `var' "`varlab' (%)"
			}
		}
		* drop uneccessary variables
		if inlist(`webopt', 2) {
			local dropvars1 = "cal1_m?_*_* cal2_m?_*_ep_cvdr cal2_m?_*_ep_cvdd" 
			local dropvars2 = "lp*_cvdd *erfc_cvdr_m? *erfc_cvdd_m?"
			local dropvars = "`dropvars1' `dropvars2'"
			drop `dropvars'
		}
		* rename variable prefixes from erfc to who
		if "`renpfix'" ~= "norenpfix" {
			renpfix erfc_ who_
			foreach i of numlist `calnos' {
				renpfix cal`i'_erfc_ cal`i'_who_
			}
		}
		gen _WHO_END = .
	}		/* end of quinoi */
	
	* tabulate recalibration parameters
	if !missing("`tabparms'") {
	*	local epvars = "ep_chdmi ep_crbv ep_cvdr ep_cvdd"
		preserve
		foreach i of numlist `calnos' {
			local cpfix = "cal`i'_"		/* 1 = country-specific, 2 = GBD regions */
			local caltype = cond(`i' == 1, "Country", "Region")
			foreach epvar in `epvars' {
				local coefvars = "`cpfix'm?_*`epvar'"
				capture ds `coefvars'
				if _rc == 0 {
					di _newline as input "`caltype' " as res "`attime'-yr risk " ///
					as input "calibration parameters for epvar = " as res "`epvar'" _c
					tabdisp calyear if !missing(_mergecal), c(`coefvars') by(ccode sex) concise
				}
			}
		}
		restore
	}
	* keep minimal variables in dataset if nodetails
	if inlist(`webopt', 1, 2) {
		qui ds _WHO_START - _WHO_END
		local dropvars = r(varlist)
		local calvars = cond(inlist(`webopt', 1), "cal?_*_*_m?", "cal?_*_cvdx_m?")
		if "`details'" ~= "details" {
			local keepvars = "_WHO_START calyear gbdreg gbdregid countryid _mergecal `calvars' _WHO_END"
			qui ds `keepvars'
			local keepvars = r(varlist)
			local dropvars: list dropvars - keepvars
			drop `dropvars'
		}
		if "`relabel'" ~= "norelabel" {
			* relabel variables with better descriptive labels
			foreach var of varlist cal?_*_*_m? {
				local varlab: variable label `var'
				local varlab: subinstr local varlab "`t'-year risk cal1" "Country `t'-yr risk ", all
				local varlab: subinstr local varlab "`t'-year risk cal2" "Region `t'-yr risk ", all
				local varlab: subinstr local varlab "CVDx" "CVD"
				local varlab: subinstr local varlab "Model 1" "lab"
				local varlab: subinstr local varlab "Model 2" "non-lab"
				local varlab: subinstr local varlab "cal1WHO-" "CouWHO-", all
				local varlab: subinstr local varlab "cal2WHO-" "RegWHO-", all
				label variable `var' "`varlab'"
			}
			if "`details'" == "details" {
				* relabel additional variables with better descriptive labels
				foreach var of varlist lp_m?_* *_*_m? {
					local varlab: variable label `var'
					if missing("`shortlab'") {
						local varlab: subinstr local varlab "Model 1" "lab"
						local varlab: subinstr local varlab "Model 2" "non-lab"
					}
					else {
						local varlab: subinstr local varlab " Model 1" "M1"
						local varlab: subinstr local varlab " Model 2" "M2"
					}
					if strmatch("`var'", "who_*_m?") {
						local varlab = cond(missing("`shortlab'"), "ERFC 10-yr risk `varlab'", "ERFC`varlab'")
					}
					label variable `var' "`varlab'"
				}
			}
		}
	}
	* describe and summarize risk variables
	des _WHO_START - _WHO_END
	preserve
	keep if !missing(_mergecal)
	bysort gbdregid: tabstat cal?_*_*_m?, by(ccode) stats(n mean sd min p25 p50 p75 max) ///
	col(stats) longstub varw(16)
	restore

end

exit		/* force exit to ignore example code below */

* examples
* CVD incidence risk calculations
whocvdrisk
whocvdrisk, attime(10)
whocvdrisk, attime(10) calyear(2017)
whocvdrisk, attime(10) shortlab
whocvdrisk, attime(10) shortlab calyear(2017) 
whocvdrisk, attime(5)
whocvdrisk, attime(5) calyear(2017)
whocvdrisk, attime(5) shortlab
whocvdrisk, attime(5) shortlab calyear(2017) 

* managing generated variables
des _WHO_START - _WHO_END
drop _WHO_START - _WHO_END
