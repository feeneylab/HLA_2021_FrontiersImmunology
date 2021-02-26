*** Location of ritest package	
	sysdir set PLUS ""

*** Set locals
	local envseed : env varSeed
	local outcome : env varOut
	local exp : env varExp
	local base : env varBase
	local no2ds2 : env varNo2DS2
	local sens : env varSens
	
***	Sensitivity analyses 
  * Sens = "" : Primary analysis
  * Sens = 1 : restrict to site==1 (Jinja only)
  * Sens = 2 : restrict to site==2 (Kanungu only)
  * Sens = 3 : restrict to site==3 (Tororo only)
  * Sens = 4 : restrict to child==1 (children only)
  * Sens = 5 : restrict to hbs_r==0 (sickle cell wildtype only)
  * Sens = 6 : para_lamp_q (LAMP included in measure of parasite prevalence)
  * Sens = 7 : HLA 2 digit exposures
  * Sens = 8 : para_lamp_q_all (Children + adults - LAMP included in measure of parasite prevalence)
  
  * Sens = *** : use site_bantu in covariate list (site/Bantu combined variable instead of siteid)
  * Sens = *** : restrict to child==0 (adults only - need to have only ind random effect - no HH random effect)
	
 ** Standard covariate list
	local cov = "male eir_ln ib3.siteid" /* Different age var in each model */	
	
 ** Set locals for if statements to run subgroups for sensitivity analyses
	if "`sens'" == "1"|"`sens'" == "2"|"`sens'" == "3" {
		local sens_note = "Sensitivity analysis: siteid==`sens'"
		local ifsite = "`sens'"
		local folder = "Sens_Site`sens'"
		local cov = "male eir_ln"
	}
	
  * These match main analysis file
	if "`sens'" == "4" {
		local sens_note = "Sensitivity analysis: Child"
		local ifchild = 1
		local folder = "Sens_Child"
	}
	
	if "`sens'" == "5" {
		local sens_note = "Sensitivity analysis: Sickle cell wildtype"
		local ifhbs = 1
		local folder = "Sens_HBS"
	}
	
	if "`sens'" == "6" {
		local sens_note = "Sensitivity analysis: LAMP parasite prevalence (Child)"
		local folder = "Sens_LAMP"
	}
	
	if "`sens'" == "7" {
		local sens_note = "Sensitivity analysis: HLA 2-digit"
		local folder = "Sens_HLA2"
	}
	
	if "`sens'" == "8" {
		local sens_note = "Sensitivity analysis: LAMP parasite prevalence (All)"
		local folder = "Sens_LAMPall"
	}
	
	/*	if "`sens'" == "***" {
		local sens_note = "Sensitivity analysis: Adult"
		local ifadult = 1
		local folder = "Sens_Adult"
	} */
	
	/* 	if "`sens'" == "***" {
		local sens_note = "Sensitivity analysis: Site_Bantu"
		local ifsitebantu = 1
		local folder = "Sens_Site_Bantu"
		local cov = "male eir_ln ib3.site_bantu"/* Different age var in each model */
	} */
	
	
 ** Site/Bantu combined variable covariate list
	if "`ifsitebantu'" == "1" {
		local cov = "male eir_ln ib3.site_bantu"/* Different age var in each model */
	}
	
*** Log/Output Settings
	capture log close
	
 ** Filename
	if "`base'"=="" {
		local filename = "`outcome'_`exp'"
	}
	
	else if "`base'"!="" {
		local filename = "`outcome'_`exp'_b`base'"
	}
	
	if "`no2ds2'"=="1" {
		
		local no2ds2_note = "No 2DS2 subset only"
		
		if "`base'"=="" {
			local filename = "`outcome'_`exp'_no2DS2"
		}
		
		else if "`base'"!="" {
			local filename = "`outcome'_`exp'_b`base'_no2DS2"
		}
	}
	
	if "`sens'"=="" {
		log using "Log/`filename'", replace
		local outcome_data	= "Output_Data/`filename'"
		local outcome_csv	= "Output/`filename'"
	}
	
	else if "`sens'"!="" {
		log using "`folder'/Log/`filename'", replace
		local outcome_data	= "`folder'/Output_Data/`filename'"
		local outcome_csv	= "`folder'/Output/`filename'"
	}
	
*** Settings
	di "Exp: `exp'"
	di "Base: `base'"
	di "Outcome: `outcome'"
	di "`no2ds2_note'"
	di "`sens_note'"
	di `envseed'
	
	local reps = 10000
	
	local seed = `envseed' + 587124
	set seed `seed'
	di `seed'
	set maxiter 500
	
*** Call Data
	local data_date = "2020-09-28"
	use "Data/Malaria_HLA_Variables_`data_date'.dta", clear		

 ** If sensitivity analysis: keep data
	if "`ifsite'" == "1"|"`ifsite'" == "2"|"`ifsite'" == "3" {
		tab siteid
		keep if siteid==`ifsite'  /* One site only */
		tab siteid
	}
 
	if "`ifchild'" == "1" {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if "`ifadult'" == "1" {
		tab child
		keep if child==0  /* Adults only */
		tab child
	}
	
	if "`ifhbs'" == "1" {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}
	
 ** Set model variable and coefficient locals
	quietly tab `exp'
	
  * Binary/Continuous
	if r(r) == 2| r(r) > 4 {
	
		local exp_m = "`exp'"
		local extract = "`outcome'_`exp' = _b[`exp']"
		
		if "`no2ds2'"=="1" {
			local extract = "`outcome'_`exp'_no2DS2 = _b[`exp']"
		}
	}
		
  * Categorical (only 3 or 4 levels)
	if r(r) == 3|r(r) == 4 {
		
		levelsof `exp', local(levels)
			
		if "`base'" == "" {
			local exp_m = "i.`exp'" 
			di "`exp_m'"
			quietly sum `exp'
			local base = r(min)
			di `base'
			}	
		
		else {
			sum `exp'
			local exp_m = "ib`base'.`exp'"
			di "`exp_m'"
			}	
			
		local levels_extract: list levels- base /* List all levels except base */
		di "`levels_extract'"

		foreach l of local levels_extract { /* Loop through all levels except base to specify what to export */
			local extract_raw "`outcome'_`exp'_b`base'_`l' = _b[`l'.`exp']"
			
			if "`no2ds2'"=="1" {
				local extract_raw "`outcome'_`exp'_b`base'_`l'_no2DS2 = _b[`l'.`exp']"
			}
			
			di "`extract_raw'"
			local extract "`extract'" "`extract_raw' "
			di "`extract'"
			}
		}
	
********************************************************************************	
	
*** Log parasite density: Malaria visits

	if "`outcome'" == "para_ln_m" {

	*** Prepare data
	 ** Keep non-routine visits
		keep if visittype == 1
		
	 ** Keep if parasite density > 0
		keep if parasitedensity>0 & parasitedensity<.
		
	 ** Drop if do NOT have malaria/do NOT have fever
		drop if febrile==0 /* 2 cases, checked with incident malaria */

*** Permutations/Model		
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mixed para_ln `exp_m' c.age##i.child `cov' || hhid: || id:, reml
		
	}
		
********************************************************************************	
		
*** Log parasite density: Routine visits	
	
	else if "`outcome'" == "para_ln_r" {
	
	*** Prepare data
	 ** Keep routine visits
		keep if visittype != 1
		
	 ** Keep if parasite density > 0
		keep if parasitedensity>0 & parasitedensity<.

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mixed para_ln `exp_m' c.age##i.child `cov' || hhid: || id:, reml		
	
	}
	
********************************************************************************		
	
*** Incident Malaria: Yearly	

	else if "`outcome'" == "inc_mal" {

	*** Prepare data
	 ** Keep annual variables
		keep id hhid site* *bantu* hla* year incidentmalaria_year pt_year age_y agecat_y male eir_ln g6pd_r hbs_r alphathal_r child
		duplicates drop
		
		drop if pt_year==0

	 ** Create offset
		gen pt_days_log = log(pt_year/365.25)

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mepoisson incidentmalaria_year `exp_m' c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
		
	}
		
********************************************************************************	

*** Parasite Prevalence: Quarter	
	
	else if "`outcome'" == "para_prev" {
	
	*** Prepare data
	 ** Keep quarterly variables
		keep id hhid site* *bantu* hla* year q para_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
		duplicates drop

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	melogit para_q `exp_m' c.age_q##i.child `cov' || hhid: || id:, or		

	}
	
********************************************************************************	

*** Probability of symptoms if infected	

	else if "`outcome'" == "psi" {

	*** Call data	
		use "Data/Malaria_HLA_Variables_PSI_`data_date'.dta", clear /* Use PSI data file */
		
	 ** If sensitivity analysis: keep data
		if "`ifsite'" == "1"|"`ifsite'" == "2"|"`ifsite'" == "3" {
			tab siteid
			keep if siteid==`ifsite'  /* One site only */
			tab siteid
		}
	 
		if "`ifchild'" == "1" {
			tab child
			keep if child==1  /* Children only */
			tab child
		}
		
		if "`ifadult'" == "1" {
			tab child
			keep if child==0  /* Adults only */
			tab child
		}
		
		if "`ifhbs'" == "1" {
			tab hbs_r
			keep if hbs_r==0 /* Wildtype sickle cell */
			tab hbs_r
		}	
		
	*** Prepare data
	 ** Data prepared in variable creation do-file

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	melogit febrile `exp_m' c.age##i.child `cov' || hhid: || id:, or
	
	}
	
********************************************************************************	
	
*** Temperature: All parasite+ visits		
	
	else if "`outcome'" == "temp" {
	
	*** Prepare data
	 ** Keep if parasite density > 0
		keep if parasitedensity>0 & parasitedensity<.

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	mixed temperature `exp_m' c.age##i.child `cov' para_ln || hhid: || id:, reml	
	
	}
	
********************************************************************************	

*** LAMP/Parasite Prevalence: Quarter (CHILDREN ONLY)	
	
	else if "`outcome'" == "lamp_prev" {
	
	*** Prepare data
	 ** Keep quarterly variables (CHILDREN ONLY)
		keep id hhid site* *bantu* hla* year q para_lamp_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
		tab child
		keep if child==1
		tab child
		duplicates drop

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	melogit para_lamp_q `exp_m' c.age_q `cov' || hhid: || id:, or		

	}
	
********************************************************************************	

*** LAMP/Parasite Prevalence: Quarter (CHILDREN + ADULTS)	
	
	else if "`outcome'" == "lamp_all" {
	
	*** Prepare data
	 ** Keep quarterly variables
		keep id hhid site* *bantu* hla* year q para_lamp_q_all age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
		tab child
		duplicates drop

*** Permutations/Model
	ritest `exp' "`extract'", reps(`reps') cluster(id) noisily reject(e(converged)!=1) saving("`outcome_data'", replace): ///
	melogit para_lamp_q_all `exp_m' c.age_q##i.child `cov' || hhid: || id:, or		

	}
	
********************************************************************************
********************************************************************************
	
*** Fill in data for export

	gen outcome = ""
	gen exp = ""
	gen coef = ""
	gen nobs = .
	gen obs = .
	gen cttrue = .
	gen repsn = .
	gen repsnm = .
	gen p = .
	gen cilb = .
	gen ciub = .
	gen seed = .
	gen base = .
	gen sens = ""
	
	local num_extract = r(k_exp)
	di `num_extract'
	
	forvalues i=1/`num_extract' { /* Loop through number of coefficients extracted */

		replace outcome = "`outcome'" if _n==`i'
		
		replace exp = "`exp'" if _n==`i'
		
		replace sens = "`sens'" if _n==`i'
		
		mat matobs = r(b)
		di "`: word `i' of `: colnames matobs''"
		replace coef = "`: word `i' of `: colnames matobs''" if _n==`i'
		
		replace nobs = r(N) if _n==`i'
		
		replace obs = matobs[1,`i'] if _n==`i'
		
		replace repsn = r(N_reps) if _n==`i'
		
		mat matct = r(c)
		replace cttrue = matct[1,`i'] if _n==`i'
		
		mat matreps = r(reps)
		replace repsnm = matreps[1,`i'] if _n==`i'
		
		mat matp = r(p)
		replace p = matp[1,`i'] if _n==`i'
		
		mat matci = r(ci)
		replace cilb = matci[1,`i'] if _n==`i'
		replace ciub = matci[2,`i'] if _n==`i'
		
		replace seed = `seed' if _n==`i'
		capture replace base = `base' if _n==`i'
		
		}
		
		keep outcome-sens
		
		keep if _n<=`num_extract'
	
*** Export data
	export delimited using "`outcome_csv'.csv", replace
	
	clear
	exit
