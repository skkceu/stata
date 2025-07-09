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
S Kaptoge: Apr 2022
* Wrapper program called by SCORE2 algorithims to display risk region information for countries.
------------------------------------------------------------------------------------------------------------*/

capture program drop score2regs
program define score2regs, rclass
	
	version 13
	syntax [, ccodes escyear(integer 2021) noWEB noEXTRA WEBVERSion(integer 2) ///
	CALRUN(integer 1) CALPARMS(string) pause]
	
	preserve
	local escdir = "\\me-filer1\groups$\MEU\ERFC\Analysis\METHODS\PubHealth\SCORE_ESC"
	local weburl = "http://ceu.phpc.cam.ac.uk/software/erfc"
	local quinoi = "qui"
	`quinoi' {
		* web settings
		local webopt = cond("`web'" == "noweb", 0, `webversion')
		local calrunparm = "CALRUN" + cond(`calrun' < 10, "0", "") + string(`calrun')
		if inlist(`webopt', 0, 1, 2) {
			if inlist(`webopt', 0) {
				local webdir = "`escdir'\WEB\score2risk\\`escyear'\\`calrunparm'"
			}
			else {
				local webdir = "`weburl'/score2risk/`escyear'/`calrunparm'"
			}
			if missing("`calparms'") {
				local calparms = "`webdir'\esc_score_recalib_parms_v05.dta"
			}
		}
		* get parameters as supplied
		use "`calparms'", clear
		tempfile calparms
		save `calparms'
		if "`pause'" == "pause" {
			pause on
			pause
		}
		* display country code mappings to GBD regions
		if inlist("`ccodes'", "ccodes", "")  {
			egen tag1 = tag(country)
			clonevar riskregid = riskreg
			label define riskregid 1 "Europe Low risk region" 2 "Europe Moderate risk region" 3 "Europe High risk region" 4 "Europe Very high risk region", modify
			label values riskregid riskregid
			label values riskreg
			noi tabdisp country if tag1 == 1, by(riskregid) concise c(ccode riskreg)
			mkmat riskreg if tag1 == 1, matrix(riskregs) rownames(country) roweq(ccode)
			return matrix riskregs = riskregs, copy
		}
		restore
		if "`extra'" == "noextra" {
			noi di as text "Done for now ..."
			exit
		}
	}
end

exit		/* force exit to ignore example code below */

* examples
* display country risk region information
score2regs
score2regs, noweb
score2regs, webversion(0)
score2regs, webversion(1)
score2regs, webversion(2)

