*** MALARIA HLA ANALYSIS
 ** Added 10-31-20 HLA class 2 homozygosity

*** SETTINGS AND PATHS
    clear all
	set more off
    capture log close
	program drop _all
	macro drop _all
	
	capture cd ""
	display _rc
	if _rc==0 {
		global jdlaptop = 1
		}
	else {
		global jdlaptop = 0
		}

	capture cd "" /* Run on MyResearch */
	display _rc

	local datapath  	= "Data/"
 	local logpath   	= "Log/"
	global outputpath 	= "Output/" /* Must be global to export in programs */
	
 ** Create date local to name files
	local currdate_raw: display %td_CCYY_NN_DD date(c(current_date), "DMY")
	local currdate = subinstr(trim("`currdate_raw'"), " " , "-", .)

	
	*local inputdate = "2020-06-25"
	local inputdate = "2020-10-31"
	
	local inputfile	 		= "Malaria_HLA_Variables_`inputdate'"
	local inputfile_psi	 	= "Malaria_HLA_Variables_PSI_`inputdate'"
	
 ** Set exposure list
	/* Set only this to 1 to run primary 4-digit HLA analyses, leave all other locals 0 */
	local explist_hla_l = 1
	
	/* Set only this to 1 to run 2-digit sensitivity HLA analyses for all outcomes, leave all other locals 0 */
	local explist_hla_s = 0
 
 ** Run models split by site (in addition to all together) if run_site = 1 for sensitivity analyses
	global run_site = 0
	
 ** Set locals for if statements to run subgroups for sensitivity analyses
	local ifchild = 0
	local ifhbs = 0
	
	/*	To run adults only:
		Random effect should be individual only, NOT household */
	
 ** Set local to change covariate list to include site_bantu instead of siteid for sensitivity analyses
	local ifsitebantu = 0
	
 ** Set to just run parasite prevalence models for sensitivity analyses
	local paraprevonly = 0
	
 ** Set to just run LAMP models for sensitivity analyses
	local lamponly = 0 /* Toggle to just run lamp models for sensitivity analyses */
	
 ** Set to just run LAMP models for sensitivity analyses
	local lamponlyall = 0 /* Toggle to just run lamp models (children and adults) for sensitivity analyses */	
	
 ** Exposures
	if `explist_hla_l' == 1 {
		local exp_list_all = "hla_l hla_homo"
	
		global excelfile	  	= "Malaria_HLA_Analysis_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_HLA_Analysis_`currdate'"
		
		local filename = "HLA"
	}
	
	if `explist_hla_s' == 1 {
		local exp_list_all = "hla_s"
	
		global excelfile	  	= "Malaria_HLA_s_Analysis_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_HLA_s_Analysis_`currdate'"
	}	
	
 ** File Names
	if $run_site == 1 {
		global excelfile	  	= "Malaria_`filename'_Analysis_allsites_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_allsites_`currdate'"
	}
	
	if `ifchild' == 1 {
		global excelfile	  	= "Malaria_`filename'_Analysis_child_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_child_`currdate'"
	}
	
	if `ifhbs' == 1 {
		global excelfile	  	= "Malaria_`filename'_Analysis_hbs_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_hbs_`currdate'"
	}
	
	if `lamponly' == 1 {
		global excelfile	  	= "Malaria_`filename'_Analysis_lamp_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_lamp_`currdate'"
	}
	
	if `lamponlyall' == 1 {
		global excelfile	  	= "Malaria_`filename'_Analysis_lampall_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_lampall_`currdate'"
	}
	
 ** Standard covariate list
	local cov = "male eir_ln ib3.siteid" /* Different age var in each model */
	
 ** Site/Bantu combined variable covariate list
	if `ifsitebantu' == 1 {
		local cov = "male eir_ln ib3.site_bantu"/* Different age var in each model */
		global excelfile	  	= "Malaria_`filename'_Analysis_sitebantu_`currdate'.xlsx"  /* Must be global to export in programs */
		local logfile	  		= "Malaria_`filename'_Analysis_sitebantu_`currdate'"
	}
	
*** OPEN LOG FILE
	log using "`logpath'/`logfile'", replace
	
 ** Set max iterations
	set maxiter 500

********************************************************************************

*** Call data
	use "`datapath'/`inputfile'", clear	

********************************************************************************	
********************************************************************************	

*** Exposures
 ** Determine which HLA alleles to model
	preserve
	keep id hla*
	duplicates drop /* Keep one copy per person to count freq of HLA to determine which to model */
	count			/* Should have 890 */

 ** HLA Locals
  * 4 Digit HLA
	local exp_hla_l = ""
	
	if `explist_hla_l' == 1 {
	foreach hla in a b c dpa1 dpb1 dqa1 dqb1 drb1 {
		foreach var of varlist hla_`hla'_l_* {
				tab `var'
				scalar `var'_tot = r(N)
				count if `var' == 1
				scalar `var'_1 = r(N)
				global `var'_p = `var'_1/`var'_tot
				di ${`var'}_p
				
				if `var'_1/`var'_tot>0.03 { /* add all HLA with > 0.03 freq to local list */
					local exp_hla_l `exp_hla_l' `var'
					di "`exp_hla_l'"
				}
			}
		}	
	}
	di "`exp_hla_l'"

  * 2 Digit HLA
	local exp_hla_s = ""
	
	if `explist_hla_s' == 1 {
	
	foreach hla in a b c dpa1 dpb1 dqa1 dqb1 drb1 {

		foreach var of varlist hla_`hla'_s_* {
				tab `var'
				scalar `var'_tot = r(N)
				count if `var' == 1
				scalar `var'_1 = r(N)
				global `var'_p = `var'_1/`var'_tot
				di ${`var'}_p
				
				if `var'_1/`var'_tot>0.03 { /* add all HLA with > 0.03 freq to local list */
					local exp_hla_s `exp_hla_s' `var'
					di "`exp_hla_s'"
				}
			}
		}	
	}
	di "`exp_hla_s'"
	
  * HLA Homo
	local exp_hla_homo ib2.hla_class1_homo_l_cat ib2.hla_class1_homo_s_cat

********************************************************************************	
********************************************************************************
	
 ** Program to create descriptive tables
	
	program desc_exp
	
		args title
		
		* Fill in table title
		putexcel A1 = "`title'"
		
		*putexcel A2 = ""
		putexcel B2 = "N"
		putexcel C2 = "%"
		putexcel D2 = "Total N"
		
		global d = 3
		
 ** Fill in table values
	foreach var of varlist $descvars {
	
		* Fill in variable label
		  local varlab: variable label `var'
		  putexcel A$d = "`varlab'"
		  
		* Total N 
		  tab `var', nol
		  scalar ntot = r(N)
		  putexcel D$d = ntot
		  
		* Count number of levels (to determine categorical or continuous)
		  local nlevels = r(r)
		  
		* Fill in for binary
		  if `nlevels'==2 {
			sum `var'
			count if `var'==r(max)
			scalar nlev = r(N)
			putexcel B$d = nlev		 /* Fill in N for that level */
			putexcel C$d = nlev/ntot*100 /* Fill in % */
			
			* Move to the next row
			  global d = $d + 1
				
			}
		  
		* Fill in for categorical {
		  else if `nlevels'>2 & `nlevels'<6 {
		  
			* Get value labels and make locals of them
			  levelsof `var', local(levels)
			  local vallab: value label `var'
			  foreach l of local levels {
				local vallab`l': label `vallab' `l'
				}
			  
			* Move to the next row
			  global d = $d + 1
		  
			* Summarize values
			  summ `var'
			  local min = r(min)
			  local max = r(max)
			  
			* Loop through each value
			  forvalues j = `min'/`max' {
		
				count if `var'==`j' /* Get N in that category */
				scalar nlev = r(N)
				putexcel B$d = nlev /* Fill in N in j category */
				putexcel C$d = nlev/ntot*100 /* Percent of total N */
				putexcel A$d ="`vallab`j''" /* Fill in value label of j */
			  
			  * Move to the next row
				global d = $d + 1
			  
			  }
			}

		}
	
	end
	
********************************************************************************	

*** Descriptive tables
	use "`datapath'/`inputfile'", clear	

	keep id hla_*_l_* hla_*_s_*
	duplicates drop
	
	if `explist_hla_l' == 1 {
		global descvars = ""
		foreach var of varlist hla_*_l_* {
			global descvars $descvars `var'
			di "$descvars"
			}
		putexcel set "${outputpath}/${excelfile}", sheet("desc_hla_l") modify
		desc_exp "HLA alleles" $descvars
	}
	
	if `explist_hla_s' == 1 {
		global descvars = ""
		foreach var of varlist hla_*_s_* {
			global descvars $descvars `var'
			di "$descvars"
			}
		putexcel set "${outputpath}/${excelfile}", sheet("desc_hla_s") modify
		desc_exp "HLA allele groups" $descvars
	}

*** Table 1
	putexcel set "${outputpath}/${excelfile}", sheet("Table 1") modify

	foreach s in 1 2 3 {
	
	local row = 2
	
	if `s'==1 {
		local col = "B"
		}
		
	if `s'==2 {
		local col = "C"
		}
		
	if `s'==3 {
		local col = "D"
		}
	
	* Household data
		use "`datapath'/`inputfile'", clear	
		keep hhid eir siteid
		duplicates drop
		
		count if siteid==`s'
		scalar n_hh = r(N)
		putexcel `col'`row' = n_hh
		local row = `row'+1
		
		sum eir if siteid==`s', d
		local median_eir = string(`r(p50)',"%9.1f")
		local min_eir = string(`r(min)',"%9.1f")
		local max_eir = string(`r(max)',"%9.1f")
		putexcel `col'`row' = "`median_eir' (`min_eir'-`max_eir')"
		local row = `row'+1
		
	* Adult data
			use "`datapath'/`inputfile'", clear	
			keep if child==0
			keep id ageyrs male siteid incidencerate
			duplicates drop
			
			count if siteid==`s'
			local n_people = r(N)
			putexcel `col'`row' = `n_people'
			local row = `row'+1
			
			gen incidencerate_yr = incidencerate*365
			sum incidencerate_yr if siteid==`s', d
			local median_ir = string(`r(p50)',"%9.1f")
			local min_ir = string(`r(min)',"%9.1f")
			local max_ir = string(`r(max)',"%9.1f")
			putexcel `col'`row' = "`median_ir' (`min_ir'-`max_ir')"
			local row = `row'+1
			
			use "`datapath'/`inputfile'", clear		
			* Keep quarterly variables
			keep if child==0
			keep id siteid para_q year q
			duplicates drop
			count if para_q==1 & siteid==`s'
			local n_para_q = r(N)
			count if siteid==`s' & para_q!=.
			local n_q = r(N)
			local para_q_p_raw = `n_para_q'/`n_q'*100
			local para_q_p = string(`para_q_p_raw',"%9.1f")
			putexcel `col'`row' = "`para_q_p'%"
			local row = `row'+1

	* Child data
			use "`datapath'/`inputfile'", clear	
			keep if child==1
			keep id ageyrs male siteid incidencerate
			duplicates drop
			
			count if siteid==`s'
			local n_people = r(N)
			putexcel `col'`row' = `n_people'
			local row = `row'+1
			
			sum ageyrs if siteid==`s', d
			local median_ageyrs = string(`r(p50)',"%9.1f")
			local min_ageyrs = string(`r(min)',"%9.1f")
			local max_ageyrs = string(`r(max)',"%9.1f")
			putexcel `col'`row' = "`median_ageyrs' (`min_ageyrs'-`max_ageyrs')"
			local row = `row'+1
			
			gen incidencerate_yr = incidencerate*365
			sum incidencerate_yr if siteid==`s', d
			local median_ir = string(`r(p50)',"%9.1f")
			local min_ir = string(`r(min)',"%9.1f")
			local max_ir = string(`r(max)',"%9.1f")
			putexcel `col'`row' = "`median_ir' (`min_ir'-`max_ir')"
			local row = `row'+1
			
			use "`datapath'/`inputfile'", clear		
			* Keep quarterly variables
			keep if child==1
			keep id siteid para_q year q
			duplicates drop
			count if para_q==1 & siteid==`s'
			local n_para_q = r(N)
			count if siteid==`s' & para_q!=.
			local n_q = r(N)
			local para_q_p_raw = `n_para_q'/`n_q'*100
			local para_q_p = string(`para_q_p_raw',"%9.1f")
			putexcel `col'`row' = "`para_q_p'%"
			local row = `row'+1		
			
		
	* All data
		use "`datapath'/`inputfile'", clear	
		keep id pt_study_month hbs hbs_r siteid
		duplicates drop
		
		sum pt_study_month if siteid==`s', d
		local median_pt = string(`r(p50)',"%9.1f")
		local min_pt = string(`r(min)',"%9.1f")
		local max_pt = string(`r(max)',"%9.1f")
		putexcel `col'`row' = "`median_pt' (`min_pt'-`max_pt')"
		local row = `row'+1
		
		local row = `row'+1
		count if siteid==`s'
		local n_tot = r(N)
		foreach i in 0 1 2 9 {
			count if hbs==`i' & siteid==`s'
			local n_hbs = r(N)
			local hbs_p_raw = `n_hbs'/`n_tot'*100
			local hbs_p = string(`hbs_p_raw',"%9.1f")
			putexcel `col'`row' = "`n_hbs' (`hbs_p'%)"
			local row = `row'+1
		}

	}
	
  * Totals for text
	use "`datapath'/`inputfile'", clear	
	keep id hhid siteid pt_study_month child
	duplicates drop
	tab child
	quietly tab hhid
	display r(r)
	sum pt_study_month
	
	use "`datapath'/`inputfile'", clear	
	keep id siteid bantu
	duplicates drop
	tab siteid bantu, row
	
  * Outcomes for text
  * Malaria incidence rate for entire study (one row per person)
	use "`datapath'/`inputfile'", clear
	keep id hhid siteid incidencerate child
	duplicates drop /* One row per person per year */
	
	gen incidencerate_yr = incidencerate*365
	
	bysort siteid: sum incidencerate_yr if child==1, d
	
	bysort siteid: sum incidencerate_yr if child==0, d
	
 ** Parasite prevalence (one row per quarter)
	use "`datapath'/`inputfile'", clear		
	
  * Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_q age_q agecat_q male eir_ln g6pd_r hbs* alphathal_r child
	duplicates drop

	bysort siteid: tab para_q if child==1
	
	bysort siteid: tab para_q if child==0

********************************************************************************

*** Program to fill in exposure across the top of each sheet
	program label_reg
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		if $i == 1 {	
			putexcel C1 = "$exp"
			}
		else if $i > 1 {
			excelcol $colnum
			local colname `r(column)'
			putexcel `colname'1 = "$exp"
		}
	}
	end

*** Program to export regression results
	program reg_exp
		args b
		scalar econverged = e(converged) /* For reg_exp_one program below */
		global outcome_var = e(depvar) /* For lincom program below - not used in this program */
		
		mat mat1 = r(table)
		mat mat2 = mat1["b",1...]\mat1["pvalue",1...]\mat1["ll",1...]\mat1["ul",1...]
		mat mat3 = mat2'
		
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify	
		if $i == 1 {
			global colnum = 1 /* Start in first column */
			excelcol $colnum /* Convert column number to letter */
			local colname `r(column)' /* Create local of column letter */
			local rownum = $rownum + 1 /* Add 1 to row to get row 3 */
			putexcel `colname'$rownum = matrix(mat3), nformat("0.000") names /* Fill in matrix with column names */
			
			* Overwrite exposure name for that model with generic name for sheet
			putexcel B`rownum' = "HLA"
			
			global colnum = $colnum + 2 /* Move two columns over from current column */
			excelcol $colnum  /* Convert column number to letter */
			local colname `r(column)' /* Create local of column letter */
			putexcel `colname'`=`rownum'-1' = "`b'" /* Fill in label for coef/HR/IRR/beta */
			local rownum = $rownum +`= rowsof(mat3)'+1 /* Add the length of the matrix + 1 to move down rows */
			putexcel B`rownum' = "N" /* Fill in labels for N, # clusters, # individuals */
			putexcel B`=`rownum'+1' = "# HH"
			putexcel B`=`rownum'+2' = "# ID"
		}
		
		else if $i > 1 {
			excelcol $colnum /* Start in current column */
			local colname `r(column)' /* Create local of column letter */
			putexcel `colname'$rownum = matrix(mat3), nformat("0.000") colnames /* Fill in matrix with ONLY column names (no row labels) */
			putexcel `colname'$rownum = "`b'" /* Fill in label for coef/HR/IRR/beta */
		}
			
		global rownum = $rownum +`= rowsof(mat3)'+1 /* Add the length of the matrix + 1 to move down rows */
		mat mat4 = [e(N)],[e(N_g)] /* Capture N and number of groups */
		mat mat5 = mat4' /* Transpose the matrix */
		putexcel `colname'$rownum = matrix(mat5) /* Fill in N and number of groups */
	}
	
	end
	
*** Program to export results in one row to merge with permutation results
	global j = 1 /* Set row counter */
	global site_id = 0 /* Set flag for site id */
	
	program reg_exp_one
	
		putexcel set "${outputpath}/${excelfile}", sheet("results_one") modify
		* Add row labels
		
		if $j == 1 {
			putexcel A1 = "outcome"
			putexcel B1 = "outcome_var"
			putexcel C1 = "exp"
			putexcel D1 = "coef"
			putexcel E1 = "p"
			putexcel F1 = "coef_ci_lb"
			putexcel G1 = "coef_ci_ub"
			putexcel H1 = "n_obs"
			putexcel I1 = "n_hh"
			putexcel J1 = "n_id"
			putexcel K1 = "site_id"
			putexcel L1 = "converged"
			putexcel M1 = "base"
			putexcel N1 = "prevalence"
			global j = $j + 1
			}
		
		di "$exp"
		
		if substr("$exp", 1, 2) == "i." {
			global exp_noi = regexr("$exp", "i\.", "")
			}
		else if substr("$exp", 1, 2) == "ib" {
			global exp_noi = regexr("$exp", "ib[0-9]\.", "")
			}
		else {
			global exp_noi = "$exp"
			}
			/* Remove i. from exposure to use string to match coefficient column names */
			
		di "$exp_noi"
		*mat mat1 = r(table) /* Use matrix from reg_exp; if don't run reg_exp first, need to unstar */
		
		local matnametest : colnames mat1 /* Save column names (coefficients) into a local */
		di "`matnametest'"
		di wordcount("`matnametest'")
		local matnametest: subinstr local matnametest "$exp_noi" "$exp_noi", all count(local c) 
			/* Use substring function to count number of coefficients with exposure name in them;
			   this equals the number of columns to extract from the matrix */
		di `c' /* Check number of columns to extract */
		
		if "$outcome_var" != "para_ln" {
			mat mat2 = mat1["b",1..`c']\mat1["pvalue",1..`c']\mat1["ll",1..`c']\mat1["ul",1..`c']
			mat mat3 = mat2'
		}
		
		if "$outcome_var" == "para_ln" {
			/* Exponeniate coefs for log parasite density for interpretability */
			mat mat2 = mat1["b",1..`c']\mat1["pvalue",1..`c']\mat1["ll",1..`c']\mat1["ul",1..`c']
			mat mat2a = mat2
			foreach i in 1 3 4 {
				/* Exponentiate rows for b, ll, ul */
					forvalues j = 1/`c' {
						mat mat2a[`i', `j'] = exp(mat2[`i', `j'])
					}
				}
			mat mat3 = mat2a'
		}
		
		putexcel A$j = "$outcome"
		putexcel B$j = matrix(mat3), rownames
			
		scalar n_obs = e(N)
		mat n = e(N_g)
		scalar n_hh = n[1,1]
		scalar n_id = n[1,2]
		
		putexcel H$j = n_obs
		putexcel I$j = n_hh
		putexcel J$j = n_id
		putexcel K$j = "$site_id"
		putexcel L$j = econverged
		if "${${exp_noi}_p}"!="" {
			di "${${exp_noi}_p}"
			putexcel N$j = ${${exp_noi}_p}
		}
		
			if `c' > 1 { /* Fill in 2nd row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				putexcel K$j = "$site_id"
				putexcel L$j = econverged
				}
			
			if `c' > 2 { /* Fill in 3rd row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				putexcel K$j = "$site_id"
				putexcel L$j = econverged
				}
				
			if `c' > 3 { /* Fill in 4th row of coefficients if more than one coefficient */
				global j = $j + 1
				putexcel A$j = "$outcome"
				putexcel H$j = n_obs
				putexcel I$j = n_hh
				putexcel J$j = n_id
				putexcel K$j = "$site_id"
				putexcel L$j = econverged
				
				}
		
			global j = $j + 1
	end
	
*** Program to export lincom results in one line
	program lin_exp_one
	
		putexcel set "${outputpath}/${excelfile}", sheet("results_one") modify
		
		if substr("${exp}",1,2)=="i." { /* For i. (all base = 0), lincom to compare level 2 to level 1 */
			global exp1 = "1." + substr("${exp}",3,.)
			global exp2 = "2." + substr("${exp}",3,.)
			scalar base = 1
			global explin = "$exp2"
			
			if "$outcome"!="temp" {
				lincom "$exp2" - "$exp1", eform
				}
			else if "$outcome"=="temp" {
				lincom "$exp2" - "$exp1"
				}
			
			putexcel A$j = "$outcome"
			putexcel B$j = "$outcome_var" /* From scalar in reg_exp program above */
			putexcel C$j = "$explin"
			
			scalar coef = r(estimate)
			putexcel D$j = coef
			scalar p_lin= r(p)
			putexcel E$j = p_lin
			scalar lb = r(lb)
			putexcel F$j = lb
			scalar ub = r(ub)
			putexcel G$j = ub

			putexcel H$j = n_obs /* Scalars for H, I, J, L from programs above */
			putexcel I$j = n_hh
			putexcel J$j = n_id
			putexcel K$j = "$site_id"
			putexcel L$j = econverged
			di econverged
			putexcel M$j = base
		
		global j = $j + 1
		}
		
		if substr("$exp", 1, 2) == "i." {
			global exp_noi = regexr("$exp", "i\.", "")
			}
		else if substr("$exp", 1, 2) == "ib" {
			global exp_noi = regexr("$exp", "ib[0-9]\.", "")
			}
		else {
			global exp_noi = "$exp"
			}
			/* Remove i. from exposure to use string to match coefficient column names */
		
		sum $exp_noi
		
		if r(max)==3 & substr("${exp}",1,2)=="i." { /* For variables with max = 3 and base = 0 */
		
			global exp3 = "3." + substr("${exp}",3,.)
			scalar base = 2
			global explin = "$exp3"
			
			if "$outcome"!="temp" {
				lincom "$exp3" - "$exp2", eform /* Lincom to compare level 3 to level 2 */
				}
			else if "$outcome"=="temp"{
				lincom "$exp3" - "$exp2" /* Lincom to compare level 3 to level 2 */
				}
		
			putexcel A$j = "$outcome"
			putexcel B$j = "$outcome_var" /* From reg_exp program above */
			putexcel C$j = "$explin"
			
			scalar coef = r(estimate)
			putexcel D$j = coef
			scalar p_lin= r(p)
			putexcel E$j = p_lin
			scalar lb = r(lb)
			putexcel F$j = lb
			scalar ub = r(ub)
			putexcel G$j = ub

			putexcel H$j = n_obs /* Scalars for H, I, J, L from programs above */
			putexcel I$j = n_hh
			putexcel J$j = n_id
			putexcel K$j = "$site_id"
			putexcel L$j = econverged
			di econverged
			putexcel M$j = base
		
		global j = $j + 1
			
			global exp3 = "3." + substr("${exp}",3,.)
			scalar base = 1
			global explin = "$exp3"
			
			if "$outcome"!="temp" {
				lincom "$exp3" - "$exp1", eform /* Lincom to compare level 3 to level 1 */
				}
			else if "$outcome"=="temp"{
				lincom "$exp3" - "$exp1" /* Lincom to compare level 3 to level 1 */
				}
				
			putexcel A$j = "$outcome"
			putexcel B$j = "$outcome_var" /* From reg_exp program above */
			putexcel C$j = "$explin"
			
			scalar coef = r(estimate)
			putexcel D$j = coef
			scalar p_lin= r(p)
			putexcel E$j = p_lin
			scalar lb = r(lb)
			putexcel F$j = lb
			scalar ub = r(ub)
			putexcel G$j = ub

			putexcel H$j = n_obs /* Scalars for H, I, J, L from programs above */
			putexcel I$j = n_hh
			putexcel J$j = n_id
			putexcel K$j = "$site_id"
			putexcel L$j = econverged
			di econverged
			putexcel M$j = base
		
		global j = $j + 1
			}
			
		if substr("${exp}",1,4)=="ib1." { /* For ib1., lincom to compare level 2 to level 0 */	
			global exp0 = "0." + substr("${exp}",5,.)
			global exp2 = "2." + substr("${exp}",5,.)
			scalar base = 0
			global explin = "$exp2"
			
			if "$outcome"!="temp" {
				lincom "$exp2" - "$exp0", eform
				}
			else if "$outcome"=="temp" {
				lincom "$exp2" - "$exp0"
				}
			
			putexcel A$j = "$outcome"
			putexcel B$j = "$outcome_var" /* From reg_exp program above */
			putexcel C$j = "$explin"
			
			scalar coef = r(estimate)
			putexcel D$j = coef
			scalar p_lin= r(p)
			putexcel E$j = p_lin
			scalar lb = r(lb)
			putexcel F$j = lb
			scalar ub = r(ub)
			putexcel G$j = ub

			putexcel H$j = n_obs /* Scalars for H, I, J, L from programs above */
			putexcel I$j = n_hh
			putexcel J$j = n_id
			putexcel K$j = "$site_id"
			putexcel L$j = econverged
			di econverged
			putexcel M$j = base
		
		global j = $j + 1
		
			}
			
		sum $exp_noi
		
		if r(max)==3 & substr("${exp}",1,4)=="ib1." { /* For ib1. with 3 levels, lincom to compare level 3 */	
			global exp3 = "3." + substr("${exp}",5,.)
			scalar base = 0
			global explin = "$exp3"
			
			if "$outcome"!="temp" {
				lincom "$exp3" - "$exp0", eform /* Lincom to compare level 3 to level 0 */
				}
			else if "$outcome"=="temp"{
				lincom "$exp3" - "$exp0" /* Lincom to compare level 3 to level 0 */
				}
		
			putexcel A$j = "$outcome"
			putexcel B$j = "$outcome_var" /* From reg_exp program above */
			putexcel C$j = "$explin"
			
			scalar coef = r(estimate)
			putexcel D$j = coef
			scalar p_lin= r(p)
			putexcel E$j = p_lin
			scalar lb = r(lb)
			putexcel F$j = lb
			scalar ub = r(ub)
			putexcel G$j = ub

			putexcel H$j = n_obs /* Scalars for H, I, J, L from programs above */
			putexcel I$j = n_hh
			putexcel J$j = n_id
			putexcel K$j = "$site_id"
			putexcel L$j = econverged
			di econverged
			putexcel M$j = base
		
		global j = $j + 1
			
			global exp3 = "3." + substr("${exp}",5,.)
			scalar base = 2
			global explin = "$exp3"
			
			if "$outcome"!="temp" {
				lincom "$exp3" - "$exp2", eform /* Lincom to compare level 3 to level 2 */
				}
			else if "$outcome"=="temp"{
				lincom "$exp3" - "$exp2" /* Lincom to compare level 3 to level 2 */
				}
				
			putexcel A$j = "$outcome"
			putexcel B$j = "$outcome_var" /* From reg_exp program above */
			putexcel C$j = "$explin"
			
			scalar coef = r(estimate)
			putexcel D$j = coef
			scalar p_lin= r(p)
			putexcel E$j = p_lin
			scalar lb = r(lb)
			putexcel F$j = lb
			scalar ub = r(ub)
			putexcel G$j = ub

			putexcel H$j = n_obs /* Scalars for H, I, J, L from programs above */
			putexcel I$j = n_hh
			putexcel J$j = n_id
			putexcel K$j = "$site_id"
			putexcel L$j = econverged
			di econverged
			putexcel M$j = base
		
		global j = $j + 1
			}
		
	end
	
*** Program to run/export AIC/BIC
	program aic_bic
	
	if $jdlaptop==1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		estat ic
		local rownum = $rownum
		putexcel B`=`rownum'+3' = "AIC"
		putexcel B`=`rownum'+4' = "BIC"
		mat mat6 = r(S)
		scalar aic = mat6[1,5]
		scalar bic = mat6[1,6]
		
		excelcol $colnum
		local colname `r(column)'
		putexcel `colname'`=`rownum'+3' = aic, nformat("0.00")
		putexcel `colname'`=`rownum'+4' = bic, nformat("0.00")
		global rownum = `rownum'+6
	}
	end
		
********************************************************************************

*** Loop through exposures
	foreach exp_list in `exp_list_all' {
		
		global explistname = "`exp_list'"
	
********************************************************************************
********************************************************************************

*** Run all outcomes if local paraprevonly = 0 and lamponly = 0 and lamponlyall = 0
	if `paraprevonly' == 0 & `lamponly' == 0 & `lamponlyall' == 0 {
	
********************************************************************************

*** Log parasite density: Routine visits
	global outcome = "para_ln_r"
	
 ** Call data
	use "`datapath'/`inputfile'", clear		

 ** Keep routine visits
	keep if visittype != 1
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.

 ** Keep if: subgroup models
	if `ifchild' == 1 {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if `ifhbs' == 1 {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}
	
 ** Empty models
	putexcel set "${outputpath}/${excelfile}", sheet(Empty) modify
	reg para_ln
	estat ic
	mixed para_ln || id:
	estat ic
	mixed para_ln || hhid: || id:
	estat ic /* Gives you AIC and BIC */

 ** Linear models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Log parasitemia: Routine visits"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
	
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		/* Crude mixed model */
		mixed para_ln `exp' || hhid: || id:, reml 
		reg_exp Coef
		aic_bic	
			
		/* Adjusted mixed model */
		mixed para_ln `exp' c.age##i.child `cov' || hhid: || id:, reml
		reg_exp Coef
		reg_exp_one
		aic_bic	
		lin_exp_one		
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				mixed para_ln `exp' c.age##i.child male eir_ln if siteid==`s' || hhid: || id:, reml
				reg_exp Coef
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			} 
		}
		
		global colnum = $colnum + 4
		global i = $i + 1	
		
	}

********************************************************************************

*** Log parasite density: Malaria care-seeking visits (non-routine)
	global outcome = "para_ln_m"
	
 ** Call data
	use "`datapath'/`inputfile'", clear		

 ** Keep non-routine visits
	keep if visittype == 1
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.
	
 ** Drop if do NOT have malaria/do NOT have fever
	drop if febrile==0 /* 2 cases, checked with incident malaria */
	
 ** Keep if: subgroup models
	if `ifchild' == 1 {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if `ifhbs' == 1 {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}	
	
 ** Empty models
	reg para_ln
	estat ic
	mixed para_ln || id:
	estat ic
	mixed para_ln || hhid: || id:
	estat ic
	
 ** Linear models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Log parasitemia: Non-routine visits"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		mixed para_ln `exp' || hhid: || id:, reml
		reg_exp Coef
		aic_bic
		
		mixed para_ln `exp' c.age##i.child `cov' || hhid: || id:, reml
		reg_exp Coef
		reg_exp_one
		aic_bic
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				mixed para_ln `exp' c.age##i.child male eir_ln if siteid==`s' || hhid: || id:, reml
				reg_exp Coef
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			} 
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
		
	}	

********************************************************************************		

*** Incident Malaria: Year
	global outcome = "inc_mal"

 ** Call data
	use "`datapath'/`inputfile'", clear	

 ** Keep annual variables
	keep id hhid site* *bantu* hla* year incidentmalaria_year pt_year age_y agecat_y male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop /* One row per person per year */
	
	drop if pt_year==0 /* Drop if person-time = 0 */

 ** Keep if: subgroup models
	if `ifchild' == 1 {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if `ifhbs' == 1 {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}
	
 ** Create offset
	gen pt_days_log = log(pt_year/365.25) /* Log person-time */
 
 ** Empty models
	poisson incidentmalaria_year, offset(pt_days_log) irr
	estat ic
	mepoisson incidentmalaria_year, offset(pt_days_log) || id:, irr
	estat ic
	mepoisson incidentmalaria_year, offset(pt_days_log) || hhid: || id:, irr
	estat ic

 ** Poisson models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Incident malaria: Annual"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		mepoisson incidentmalaria_year `exp', offset(pt_days_log) || hhid: || id:, irr
		reg_exp IRR
		aic_bic
		
		mepoisson incidentmalaria_year `exp' c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
		reg_exp IRR
		reg_exp_one
		aic_bic
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				mepoisson incidentmalaria_year `exp' c.age_y##i.child male eir_ln if siteid==`s', offset(pt_days_log) || hhid: || id:, irr
				reg_exp IRR
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			} 
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
	
	}

********************************************************************************	

*** Probability of symptoms if infected
	global outcome = "psi"

 ** Call data
	use "`datapath'/`inputfile_psi'", clear		
	
 ** Keep if: subgroup models
	if `ifchild' == 1 {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if `ifhbs' == 1 {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}
	
 ** Empty models	
	logistic febrile
	estat ic
	melogit febrile || id:, or
	estat ic
	melogit febrile || hhid: || id:, or
	estat ic
	
 ** Logistic models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Probability of symptoms if infected"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		melogit febrile `exp' || hhid: || id:, or
		reg_exp OR
		aic_bic
		
		melogit febrile `exp' c.age##i.child `cov' || hhid: || id:, or
		reg_exp OR
		reg_exp_one
		aic_bic
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				melogit febrile `exp' c.age##i.child male eir_ln if siteid==`s' || hhid: || id:, or
				reg_exp OR
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			} 
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
		
	}

********************************************************************************

*** Temperature: All parasite+ visits
	global outcome = "temp"

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.

 ** Keep if: subgroup models
	if `ifchild' == 1 {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if `ifhbs' == 1 {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}	
	
 ** Empty models
	putexcel set "${outputpath}/${excelfile}", sheet(Empty) modify
	reg temperature
	estat ic
	mixed temperature || id:
	estat ic
	mixed temperature || hhid: || id:
	estat ic

 ** Linear models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Temperature: All parasite+ visits"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
	
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		mixed temperature `exp' para_ln || hhid: || id:, reml
		reg_exp Coef
		aic_bic	
		
		mixed temperature `exp' c.age##i.child `cov' para_ln || hhid: || id:, reml
		reg_exp Coef
		reg_exp_one
		aic_bic	
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				mixed temperature `exp' c.age##i.child male eir_ln para_ln if siteid==`s' || hhid: || id:, reml
				reg_exp Coef
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			} 
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
		
	}	

} /* Close loop for running all outcomes */

********************************************************************************
********************************************************************************

*** Run LAMP if lamponly = 1
	if `lamponly' == 1	{

*** Parasite Prevalence/LAMP: Quarter (CHILDREN ONLY)
	global outcome = "lamp_prev"

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_lamp_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	keep if child==1
	duplicates drop
	
 ** Empty models
	logistic para_lamp_q
	estat ic
	melogit para_lamp_q || id:, or
	estat ic
	melogit para_lamp_q || hhid: || id:, or
	estat ic
	
 ** Logistic models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "LAMP parasite prevalence: Quarterly"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		melogit para_lamp_q `exp' || hhid: || id:, or
		reg_exp OR
		aic_bic
		
		melogit para_lamp_q `exp' c.age_q `cov' || hhid: || id:, or
		reg_exp OR
		reg_exp_one
		aic_bic
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				melogit para_lamp_q `exp' c.age_q male eir_ln if siteid==`s' || hhid: || id:, or
				reg_exp OR
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			}
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
	}	
} /* Close loop for LAMP CHILDREN models */

********************************************************************************

*** Run LAMP (ALL) if lamponlyall = 1
	if `lamponlyall' == 1	{
	
*** Parasite Prevalence/LAMP: Quarter (CHILDREN + ADULTS)
	
	global outcome = "lamp_all"

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_lamp_q_all age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop
	
 ** Empty models
	logistic para_lamp_q_all
	estat ic
	melogit para_lamp_q_all || id:, or
	estat ic
	melogit para_lamp_q_all || hhid: || id:, or
	estat ic
	
 ** Logistic models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "LAMP parasite prevalence: Quarterly"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		melogit para_lamp_q_all `exp' || hhid: || id:, or
		reg_exp OR
		aic_bic
		
		melogit para_lamp_q_all `exp' c.age_q##i.child `cov' || hhid: || id:, or
		reg_exp OR
		reg_exp_one
		aic_bic
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				melogit para_lamp_q_all `exp' c.age_q##i.child male eir_ln if siteid==`s' || hhid: || id:, or
				reg_exp OR
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			}
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
		
	}	
	
} /* Close loop for LAMP ALL models */

********************************************************************************
********************************************************************************	

*** Run if lamponly = 0 and lamponlyall = 0 (run for all parasite prevalence models)

	if `lamponly' == 0 & `lamponlyall' == 0	{
	
*** Parasite Prevalence: Quarter
	global outcome = "para_prev"

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop /* One row per person per quarter */

 ** Keep if: subgroup models
	if `ifchild' == 1 {
		tab child
		keep if child==1  /* Children only */
		tab child
	}
	
	if `ifhbs' == 1 {
		tab hbs_r
		keep if hbs_r==0 /* Wildtype sickle cell */
		tab hbs_r
	}	
	
 ** Empty models
	logistic para_q
	estat ic
	melogit para_q || id:, or
	estat ic
	melogit para_q || hhid: || id:, or
	estat ic
	
 ** Logistic models
	if $jdlaptop == 1 {
		putexcel set "${outputpath}/${excelfile}", sheet("${outcome}_${explistname}") modify
		putexcel A1 = "Parasite prevalence: Quarterly"
	}
	
	global colnum = 1
	global i = 1
	
	foreach exp in `exp_`exp_list'' {
		
		global rownum = 2
		global exp = "`exp'"
		label_reg
		
		melogit para_q `exp' || hhid: || id:, or
		reg_exp OR
		aic_bic
		
		melogit para_q `exp' c.age_q##i.child `cov' || hhid: || id:, or
		reg_exp OR
		reg_exp_one
		aic_bic
		lin_exp_one
		
		if $run_site == 1 {
			forvalues s = 1/3 {
				global site_id = `s'
				melogit para_q `exp' c.age_q##i.child male eir_ln if siteid==`s' || hhid: || id:, or
				reg_exp OR
				reg_exp_one
				aic_bic	
				lin_exp_one
				
				global site_id = 0
			}
		}
		
		global colnum = $colnum + 4
		global i = $i + 1
		
	}
}
	
*** Close loop through exposure lists	
}

********************************************************************************
********************************************************************************
********************************************************************************	

*** Trend tests for homozygosity

*** Log parasite density: Routine visits

 ** Call data
	use "`datapath'/`inputfile'", clear		

 ** Keep routine visits
	keep if visittype != 1
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.

	* Class I Trend	
	mixed para_ln ib2.hla_class1_homo_l_cat c.age##i.child `cov' || hhid: || id:, reml
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary
	mixed para_ln hla_class1_homo_l_b c.age##i.child `cov' || hhid: || id:, reml
	
	* Class II Trend
	mixed para_ln ib3.hla_class2_homo_l_cat c.age##i.child `cov' || hhid: || id:, reml
	contrast q.hla_class2_homo_l_cat, noeffects
		
	* Class II Binary
	mixed para_ln hla_class2_homo_l_b c.age##i.child `cov' || hhid: || id:, reml
		
********************************************************************************

*** Log parasite density: Malaria care-seeking visits (non-routine)
	
 ** Call data
	use "`datapath'/`inputfile'", clear		

 ** Keep non-routine visits
	keep if visittype == 1
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.
	
 ** Drop if do NOT have malaria/do NOT have fever
	drop if febrile==0 /* 2 cases, checked with incident malaria */
	
	* Class I Trend		
	mixed para_ln ib2.hla_class1_homo_l_cat c.age##i.child `cov' || hhid: || id:, reml
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary
	mixed para_ln hla_class1_homo_l_b c.age##i.child `cov' || hhid: || id:, reml
	
	* Class II Trend
	mixed para_ln ib3.hla_class2_homo_l_cat c.age##i.child `cov' || hhid: || id:, reml
	contrast q.hla_class2_homo_l_cat, noeffects
	
	* Class II Binary
	mixed para_ln hla_class2_homo_l_b c.age##i.child `cov' || hhid: || id:, reml
	
********************************************************************************		

*** Incident Malaria: Year

 ** Call data
	use "`datapath'/`inputfile'", clear	

 ** Keep annual variables
	keep id hhid site* *bantu* hla* year incidentmalaria_year pt_year age_y agecat_y male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop /* One row per person per year */
	
	drop if pt_year==0 /* Drop if person-time = 0 */

 ** Create offset
	gen pt_days_log = log(pt_year/365.25) /* Log person-time */
 
	* Class I Trend	 
	mepoisson incidentmalaria_year ib2.hla_class1_homo_l_cat c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary 
	mepoisson incidentmalaria_year hla_class1_homo_l_b c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
	
	* Class II Trend	 
	mepoisson incidentmalaria_year ib3.hla_class2_homo_l_cat c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
	contrast q.hla_class2_homo_l_cat, noeffects
	
	* Class II Binary	 
	mepoisson incidentmalaria_year hla_class2_homo_l_b c.age_y##i.child `cov', offset(pt_days_log) || hhid: || id:, irr
	

********************************************************************************	

*** Probability of symptoms if infected

 ** Call data
	use "`datapath'/`inputfile_psi'", clear		
	
	* Class I Trend	 
	melogit febrile ib2.hla_class1_homo_l_cat c.age##i.child `cov' || hhid: || id:, or
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary 
	melogit febrile hla_class1_homo_l_b c.age##i.child `cov' || hhid: || id:, or

	* Class II Trend	 
	melogit febrile ib3.hla_class2_homo_l_cat c.age##i.child `cov' || hhid: || id:, or
	contrast q.hla_class2_homo_l_cat, noeffects
	
	* Class II Binary 
	melogit febrile hla_class2_homo_l_b c.age##i.child `cov' || hhid: || id:, or

********************************************************************************

*** Temperature: All parasite+ visits
	
 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep if parasite density > 0
	keep if parasitedensity>0 & parasitedensity<.

	* Class I Trend	
	mixed temperature ib2.hla_class1_homo_l_cat c.age##i.child `cov' para_ln || hhid: || id:, reml	
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary 
	mixed temperature hla_class1_homo_l_b c.age##i.child `cov' para_ln || hhid: || id:, reml	

	* Class II Trend	
	mixed temperature ib3.hla_class2_homo_l_cat c.age##i.child `cov' para_ln || hhid: || id:, reml	
	contrast q.hla_class2_homo_l_cat, noeffects
	
	* Class II Binary 
	mixed temperature hla_class2_homo_l_b c.age##i.child `cov' para_ln || hhid: || id:, reml	

********************************************************************************

*** Parasite Prevalence/LAMP: Quarter (CHILDREN ONLY)

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_lamp_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	keep if child==1
	duplicates drop
	
	* Class I Trend	
	melogit para_lamp_q ib2.hla_class1_homo_l_cat c.age_q `cov' || hhid: || id:, or
	contrast q.hla_class1_homo_l_cat, noeffects

	* Class I Binary 
	melogit para_lamp_q hla_class1_homo_l_b c.age_q `cov' || hhid: || id:, or
	
	* Class II Trend	
	melogit para_lamp_q ib3.hla_class2_homo_l_cat c.age_q `cov' || hhid: || id:, or
	contrast q.hla_class2_homo_l_cat, noeffects

	* Class II Binary 
	melogit para_lamp_q hla_class2_homo_l_b c.age_q `cov' || hhid: || id:, or	
	
********************************************************************************

*** Parasite Prevalence/LAMP: Quarter (CHILDREN + ADULTS)
	
 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_lamp_q_all age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop
	
	* Class I Trend	
	melogit para_lamp_q_all ib2.hla_class1_homo_l_cat c.age_q##i.child `cov' || hhid: || id:, or
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary 
	melogit para_lamp_q_all hla_class1_homo_l_b c.age_q##i.child `cov' || hhid: || id:, or
	
	* Class II Trend	
	melogit para_lamp_q_all ib3.hla_class2_homo_l_cat c.age_q##i.child `cov' || hhid: || id:, or
	contrast q.hla_class2_homo_l_cat, noeffects
	
	* Class II Binary 
	melogit para_lamp_q_all hla_class2_homo_l_b c.age_q##i.child `cov' || hhid: || id:, or
	
********************************************************************************	

*** Parasite Prevalence: Quarter

 ** Call data
	use "`datapath'/`inputfile'", clear		
	
 ** Keep quarterly variables
	keep id hhid site* *bantu* hla* year q para_q age_q agecat_q male eir_ln g6pd_r hbs_r alphathal_r child
	duplicates drop /* One row per person per quarter */

	* Class I Trend		
	melogit para_q ib2.hla_class1_homo_l_cat c.age_q##i.child `cov' || hhid: || id:, or
	contrast q.hla_class1_homo_l_cat, noeffects
	
	* Class I Binary	
	melogit para_q hla_class1_homo_l_b c.age_q##i.child `cov' || hhid: || id:, or

	* Class II Trend		
	melogit para_q ib3.hla_class2_homo_l_cat c.age_q##i.child `cov' || hhid: || id:, or
	contrast q.hla_class2_homo_l_cat, noeffects
	
	* Class II Binary	
	melogit para_q hla_class2_homo_l_b c.age_q##i.child `cov' || hhid: || id:, or
	
********************************************************************************
********************************************************************************	

*** Person-time for adult LAMP vs. total	
	use "`datapath'/`inputfile'", clear		
	
	keep id pt_study_month* siteid
	duplicates drop
	
	sum pt_study_month if siteid==1, d
	sum pt_study_month_lampall if siteid==1, d
	
	sum pt_study_month if siteid==2, d
	sum pt_study_month_lampall if siteid==2, d
	
	sum pt_study_month if siteid==3, d
	sum pt_study_month_lampall if siteid==3, d
	
