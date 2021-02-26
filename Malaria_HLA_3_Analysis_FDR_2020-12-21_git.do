*** MALARIA KIR/HLA ANALYSIS

*** SETTINGS AND PATHS
    clear
	clear matrix
	clear mata
	set more off
	set maxvar 20000
    capture log close
	program drop _all
	
	local path			= ""
	
	global datapath  	= "`path'/Data"
 	local logpath   	= "`path'/Log"
	global outputpath 	= "`path'/Output"
	
 ** Create date local to name files
	local currdate_raw: display %td_CCYY_NN_DD date(c(current_date), "DMY")
	local currdate = subinstr(trim("`currdate_raw'"), " " , "-", .)

	global outputfile_all	= "Malaria_HLA_Analysis_FDR_`currdate'"
	
	local logfile	  		= "Malaria_HLA_Analysis_FDR_`currdate'"

*** OPEN LOG FILE
	log using "`logpath'/`logfile'", replace
	
********************************************************************************

*** Create program to merge original and permuted models

program merge_o_p 
	args inputfile_op inputfile_bp outputfile outcomes sheetname sites short
	
*** Call data: Original model output
	import excel "${outputpath}/`inputfile_op'", clear sheet("results_one") firstrow

 ** Keep only needed sites for sensitvity analysis by site
	if `sites' != 0 {
		keep if site_id=="`sites'"
		}
	
	gen coef_name = outcome + "_" + exp
	
	destring site_id, replace

 ** Save data 
	save "${datapath}/`outputfile'", replace

********************************************************************************

*** Call data: Batch permutation output
	import delimited "${outputpath}/`inputfile_bp'", clear varnames(1)

 ** Rename
	rename * bp_*
	rename bp_outcome outcome
	rename bp_exp exp
	rename bp_coef coef_name

 ** Merge to original model output
	merge 1:1 coef_name using "${datapath}/`outputfile'"
	drop if regexm(exp, "_homo_")==1 /* HLA homozygosity should be dropped */
	
	count if _m==1|_m==2
	local mergeerror = r(N)
	if `mergeerror'>0 {
		di in red "ERROR: `mergeerror' DID NOT MERGE FULLY"
		list coef_name if _m != 3
	}

 ** FDR of original p
  * Split p-values by outcome
  * Adjust using simes method (Simes [1986]; Benjamini and Hochberg [1995]; Benjamini and Yekutieli[2001, first method])
  	
	gen p_q = . /* Create empty variable for q-values to be combined into */
	
	foreach x in `outcomes' {
		gen p_`x' = p if outcome=="`x'"
		qqvalue p_`x', method(simes) qvalue(p_q_`x')
		count if p_q_`x'<0.05
		replace p_q = p_q_`x' if outcome=="`x'" /* Combine into one variable */
	}
	
	*browse p bp_p p_* p_q
 	
	gen bp_q = . /* Create empty variable for q-values to be combined into */
	
	foreach x in `outcomes' {
		gen bp_p_`x' = bp_p if outcome=="`x'"
		replace bp_p_`x'=1/bp_repsnm if bp_p_`x'==0
		/* Compute q value assuming bp_p = 1/non-missing if bp_p==0, most conservative assumption */
		qqvalue bp_p_`x', method(simes) qvalue(bp_q_`x')
		count if bp_q_`x'<0.05	
		replace bp_q = bp_q_`x' if outcome=="`x'" /* Combine into one variable */
	}
	
	*browse coef_name outcome p bp_p bp_q bp_*

 ** Label for export to FDR unformatted file
	* For long alleles
	if `short' == 0 {
		gen exp_r=subinstr(exp, "_l_", "*", .)
		replace exp_r=upper(subinstr(exp_r, "hla_", "", .))
		}
	
	* For short alleles
	if `short' == 1 {
		gen exp_r=subinstr(exp, "_s_", "*", .)
		replace exp_r=upper(subinstr(exp_r, "hla_", "", .))
		}
	
	* Label variables
	label var exp_r "HLA allele"
	label var p "P-value"
	label var bp_p "Empirical p-value"
	label var bp_q "FDR q-value"
	label var prevalence "Prevalence"	
	
 ** Save data 
	sort outcome exp

	export excel outcome exp_r coef p bp_p bp_q prevalence n_id using "${outputpath}/${outputfile_all}.xlsx", sheetreplace sheet(`sheetname') firstrow(variables)
end 

********************************************************************************
********************************************************************************

*** Program to fill in tables
 ** Fills in coefficient, p, bp_p, q
	program fill_num

		* Coefficient
		sum coef if coef_name == "$coefname"
		local mean = r(mean)
		if r(N)==1 {
			putexcel B$row  = `mean'
			}
		else {
			di in red "ERROR $coefname" 
			exit
			}
		* P-value
		sum p if coef_name == "$coefname"
		local mean = r(mean)
		if r(N)==1 {
			putexcel C$row  = `mean'
			}
		else {
			di in red "ERROR $coefname" 
			exit
			}
		* Permuted P-value
		sum bp_p if coef_name == "$coefname"
		local mean = r(mean)
		if r(N)==1 {
			putexcel D$row  = `mean'
			}
		else {
			di in red "ERROR $coefname" 
			exit
			}
		* FDR q-value
		sum bp_q if coef_name == "$coefname"
		local mean = r(mean)
		if r(N)==1 {
			putexcel E$row  = `mean'
			}
		else {
			di in red "ERROR $coefname" 
			exit
			}
		
		global row = $row + 1 
	end 
	
*** Fill in N for single level categorical variable
	program fill_n
		* N
		sum n_id if coef_name == "$coefname"
		local mean = r(mean)
		if r(N)==1 {
			putexcel F$row  = `mean'
			}
		else {
			di in red "ERROR $coefname" 
			exit
			}
	end
	
*** Fill in N for multiple level categorical variable
	program fill_n_mult
		sum n_id if exp == "$exp" & outcome == "$outcome"
		local mean = r(mean)
		if r(min) == r(max) {
			putexcel F$row  = `mean'
			}
		else {
			di in red "ERROR $exp $outcome"
			exit
			}
		global row = $row + 1
	end

********************************************************************************	
	
*** Tables for primary analysis
	
********************************************************************************
	
*** Program to fill out entire table
	
program fill_hla	
	args short
	putexcel set "${outputpath}/${outputtable}", sheet("Table $tablenum") modify
	putexcel A1 = "Table ${tablenum}. $tablename"
	putexcel A2 = "HLA Allele"
		if `short' == 1 {
		putexcel A2 = "HLA Allele Group"
		}
	putexcel B2 = "$coeftype"
	putexcel C2 = "P-value"
	putexcel D2 = "Empirical p-value"
	putexcel E2 = "Empirical FDR q-value"
	putexcel F2 = "N"
	
	global row = 3
	
 ** HLA
	levelsof exp, local(levels) 
	
	foreach l of local levels {
		local exp_r=subinstr("`l'", "_l_", "*", .)
		local exp_r=subinstr("`exp_r'", "_s_", "*", .)
		local exp_r=upper(subinstr("`exp_r'", "hla_", "", .))
		global coefname = "$outcome" + "_" + "`l'"
		* Exposure
		putexcel A$row  = "`exp_r'"
		* Program to fill in coef, p, bp p, q
		fill_n
		fill_num
	}
	
end

********************************************************************************
********************************************************************************

*** HLA alleles (4 digit primary analysis)

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_2020-07-07.xlsx" ///
	/*inputfile_bp*/	"hla_output_2020-07-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`currdate'" ///
	/*outcomes*/ 		"inc_mal para_ln_m para_ln_r para_prev psi temp febrile" ///
	/*sheetname*/ 		"hla_l" ///
	/*sites*/			0 ///
	/*short*/			0
	
	global outputtable	= "Malaria_HLA_SuppTables_`currdate'"
	
*** Parasite Prevalence

	global outcome = "para_prev"
	
	global tablenum = "2pp"
	global tablename = "Parasite prevalence"
	global coeftype = "Odds ratio"
	
	fill_hla 0	

*** Incident Malaria

	global outcome = "inc_mal"
	
	global tablenum = "2incmal"
	global tablename = "Incident malaria"
	global coeftype = "Incident rate ratio"
	
	fill_hla 0	
	
*** Parasite Density - Routine Visits

	global outcome = "para_ln_r"
	
	global tablenum = "2b"
	global tablename = "Parasite density: Routine quarterly parasite-positive visits"
	global coeftype = "Fold increase in parasite density"
	
	fill_hla 0	
	
*** Parasite Density - Malaria Visits

	global outcome = "para_ln_m"
	
	global tablenum = "2c"
	global tablename = "Parasite density: Malaria visits"
	global coeftype = "Fold increase in parasite density"
	
	fill_hla 0	
	
*** Probability of Symptoms if Infected

	global outcome = "psi"
	
	global tablenum = "2a"
	global tablename = "Probability of symptoms if infected"
	global coeftype = "Odds ratio"
	
	fill_hla 0	
	
*** Temperature

	global outcome = "temp"
	
	global tablenum = "2d"
	global tablename = "Temperature at parasitemic visits"
	global coeftype = "Coefficient"
	
	fill_hla 0	

********************************************************************************
********************************************************************************

*** Other Sensitivity Analyses

*** LAMP (All)	
	local sens = "lampall"

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_`sens'_2020-09-28.xlsx" ///
	/*inputfile_bp*/	"hla_output_`sens'_2020-10-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`sens'_`currdate'" ///
	/*outcomes*/ 		"lamp_all" ///
	/*sheetname*/ 		"`sens'" ///
	/*sites*/ 			0 ///
	/*short*/			0
	
		/* These files have changed LAMP all date range to include all child LAMP 
		data from entire study, with adult data only for applicable date range */

*** Parasite Prevalence

	global outcome = "lamp_all"
	
	global tablenum = "3"
	global tablename = "Parasite prevalence defined by microscopy and LAMP among children and adults"
	global coeftype = "Odds ratio"
	
	fill_hla 0	
	
********************************************************************************

*** Site 1	
	local sens = "site1"

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_allsites_2020-07-07.xlsx" ///
	/*inputfile_bp*/	"hla_output_`sens'_2020-07-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`sens'_`currdate'" ///
	/*outcomes*/ 		"para_prev" ///
	/*sheetname*/ 		"`sens'" ///
	/*sites*/ 			1 ///
	/*short*/			0

*** Parasite Prevalence

	global outcome = "para_prev"
	
	global tablenum = "4a"
	global tablename = "Parasite prevalence among residents of Jinja"
	global coeftype = "Odds ratio"
	
	fill_hla 0
	
	/* 	NOTE: hla_dpb1_l_3001 HAS TO BE MANUALLY FIXED IN EXCEL - UNABLE TO RUN
		IN SITE 1 DUE TO ONLY 1 PERSON WITH ALLELE */

********************************************************************************

*** Site 2	
	local sens = "site2"

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_allsites_2020-07-07.xlsx" ///
	/*inputfile_bp*/	"hla_output_`sens'_2020-07-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`sens'_`currdate'" ///
	/*outcomes*/ 		"para_prev" ///
	/*sheetname*/ 		"`sens'" ///
	/*sites*/ 			2 ///
	/*short*/			0

*** Parasite Prevalence

	global outcome = "para_prev"
	
	global tablenum = "4b"
	global tablename = "Parasite prevalence among residents of Kanungu"
	global coeftype = "Odds ratio"
	
	fill_hla 0	
	
********************************************************************************

*** Site 3
	local sens = "site3"

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_allsites_2020-07-07.xlsx" ///
	/*inputfile_bp*/	"hla_output_`sens'_2020-07-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`sens'_`currdate'" ///
	/*outcomes*/ 		"para_prev" ///
	/*sheetname*/ 		"`sens'" ///
	/*sites*/ 			3 ///
	/*short*/			0

*** Parasite Prevalence

	global outcome = "para_prev"
	
	global tablenum = "4c"
	global tablename = "Parasite prevalence among residents of Tororo"
	global coeftype = "Odds ratio"
	
	fill_hla 0	
	
********************************************************************************		

*** Child	

	local sens = "child"

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_`sens'_2020-07-07.xlsx" ///
	/*inputfile_bp*/	"hla_output_`sens'_2020-07-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`sens'_`currdate'" ///
	/*outcomes*/ 		"para_prev" ///
	/*sheetname*/ 		"child" ///
	/*sites*/ 			0 ///
	/*short*/			0

*** Parasite Prevalence

	global outcome = "para_prev"
	
	global tablenum = "5"
	global tablename = "Parasite prevalence among children"
	global coeftype = "Odds ratio"
	
	fill_hla 0	
	
********************************************************************************	

*** Hbs	
	local sens = "hbs"

*** Merge data
	merge_o_p ///
	/*inputfile_op*/	"Malaria_HLA_Analysis_`sens'_2020-07-07.xlsx" ///
	/*inputfile_bp*/	"hla_output_`sens'_2020-07-07.csv" ///
	/*outputfile*/ 		"Malaria_HLA_Analysis_FDR_`sens'_`currdate'" ///
	/*outcomes*/ 		"para_prev" ///
	/*sheetname*/ 		"`sens'" ///
	/*sites*/ 			0 ///
	/*short*/			0

*** Parasite Prevalence

	global outcome = "para_prev"
	
	global tablenum = "6"
	global tablename = "Parasite prevalence among HbS wildtype"
	global coeftype = "Odds ratio"
	
	fill_hla 0	
	

