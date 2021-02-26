*** MALARIA KIR/HLA VARIABLE CREATION
 ** Added 10-31-20: New dichotomous homozygosity variable and class II homozygosity 

*** SETTINGS
    clear
	clear matrix
	clear mata
	set more off
	set maxvar 20000
    capture log close
	
 ** Create date local to name files
	local currdate_raw: display %td_CCYY_NN_DD date(c(current_date), "DMY")
	local currdate = subinstr(trim("`currdate_raw'"), " " , "-", .)
	
*** PATHS
	local path			= ""
	
	local datapath  	= "`path'/Data"
 	local logpath   	= "`path'/Log"
	local outputpath 	= "`path'/Output"
	
    local inputfile_HLA			= "#694 HLA typing results 7-18-17"
	local inputfile_KIR			= "feeney_malaria KIR"
	local inputfile_c_exp		= "HLACfourdigitimputedexpression"
	local inputfile_a_exp		= "AexprJune26_VN"
	local inputfile_EIR			= "CDC light trap data all 3 sites through June 2016 FINAL"
	local inputfile_ind			= "PRISM cohort study individual study participant database through June 30th 2016 FINAL"
	local inputfile_clin 		= "PRISM cohort study all clinic visits database through June 30th 2016 FINAL"
	local inputfile_langJIN		= "Jinja PRISM1 cohort - Primary languages.xls"
	local inputfile_langKAN		= "Kanungu PRISM1-primary languages.xlsx"
	local inputfile_langTOR 	= "Tororo PRISM1 cohort - Primary languages.xls"
	
	local outputfile	 		= "Malaria_HLA_Variables_`currdate'"
	local outputfile_psi	 	= "Malaria_HLA_Variables_PSI_`currdate'"
	
	local logfile	  			= "Malaria_HLA_Variables_`currdate'"

*** OPEN LOG FILE
	log using "`logpath'/`logfile'", replace
	
********************************************************************************

*** HLA
	import excel "Data/`inputfile_HLA'.xlsx", sheet("Sheet1") firstrow case(lower) clear
	
	foreach var of varlist a_1-drb5_2 {
		gen `var's = substr(`var',1,2)
		destring `var's, replace
		
		gen hla_`var'_n = regexm(`var', "N") /* Tag if null allele */
		replace `var's=. if hla_`var'_n==1 /* Make missing if null allele */
		
		rename `var's hla_`var'_s
		gen `var'_length = length(`var')
		gen hla_`var'_l = `var'
		}

	drop if pid == "MF3273-4YNA"|pid == "MF3273-5A7P" /* UNSURE WHICH ARE THE TRUE RESULTS, DROP PER PERRI */
	drop if pid == "MF3287" /* NOT IN CLINICAL DATASET */

********************************************************************************
********************************************************************************

*** CREATE VARIABLES FOR HLA ALLELES
 ** TWO DIGITS
	foreach h in hla_a hla_b hla_c hla_dpa1 hla_dpb1 hla_dqa1 hla_dqb1 hla_drb1 {
		levelsof `h'_1_s, local(levels) 		/* get levels of HLA */
		levelsof `h'_2_s, local(levels2)
		
		foreach l of local levels {
			gen `h'_s_`l'=.						/* create empty variable */
			}
		foreach l of local levels2 {
			capture confirm variable `h'_s_`l'
			if _rc {
				gen `h'_s_`l'=.					/* create var if not already present */
				local levels `levels' `l'		/* update local with levels missing from 1 */
			}
		}
		
		foreach l of local levels {
			replace `h'_s_`l' = 1 if `h'_1_s==`l'|`h'_2_s==`l'
			* Make 0 if hla_s not missing or there is an N in the allele (null allele not expressed)
			replace `h'_s_`l' = 0 if 	`h'_1_s!=`l' & (`h'_1_s<.|`h'_1_n==1) & ///
										`h'_2_s!=`l' & (`h'_2_s<.|`h'_2_n==1)
			local hla_varlab = upper(substr("`h'",5,.)) + "*" + "`l'"
			label var `h'_s_`l' "`hla_varlab'"
			}
		browse `h'_*
	}

 ** FOUR DIGITS
	foreach h in hla_a hla_b hla_c hla_dpa1 hla_dpb1 hla_dqa1 hla_dqb1 hla_drb1 {
		local hla_length_1 = substr("`h'", 5, .) + "_1_length"
		di "`hla_length_1'"
		
		local hla_length_2 = substr("`h'", 5, .) + "_2_length"
		di "`hla_length_2'"
		
		gen `h'_1_la=subinstr(`h'_1_l, ":", "", .) if `hla_length_1' > 2 /* remove : */
		gen `h'_2_la=subinstr(`h'_2_l, ":", "", .) if `hla_length_2' > 2 /* Send long version to missing if only 2 digits */
	
		levelsof `h'_1_la, local(levels) 			/* get levels of HLA */
		levelsof `h'_2_la, local(levels2)
		
		foreach l of local levels {
			gen `h'_l_`l'=.							/* create empty variable */
			}
		foreach l of local levels2 {
			capture confirm variable `h'_l_`l'
			if _rc {
				gen `h'_l_`l'=.						/* create var if not already present */
				local levels "`levels' `l'"			/* update local with levels missing from 1 */
			}
		}
		
		foreach l of local levels {
			* Make 1 if 4-digit allele 1 or 2 matches
			replace `h'_l_`l' = 1 if `h'_1_la=="`l'"|`h'_2_la=="`l'"
			* Make 0 if 4-digit allele 1 and 2 do not match and are not missing
			replace `h'_l_`l' = 0 if `h'_1_la!="`l'" & `h'_1_la!="" & `h'_2_la!="`l'" & `h'_2_la!=""
			* Make 0 if missing 4-digit allele, but non-missing 2-digit allele does not match
			replace `h'_l_`l' = 0 if `h'_l_`l'== . & ///
					(substr("`l'",1,2) != substr(`h'_1_l,1,2)|(`hla_length_1' > 2 & `h'_1_la!="`l'")) & `h'_1_l!= "" & ///
					(substr("`l'",1,2) != substr(`h'_2_l,1,2)|(`hla_length_2' > 2 & `h'_2_la!="`l'")) & `h'_2_l!= ""
			local hla_varlab = upper(substr("`h'",5,.)) + "*" + substr("`l'",1,2) + ":" + substr("`l'",3,.)
			label var `h'_l_`l' "`hla_varlab'"
			}
		browse `h'_*
	}

********************************************************************************

 ** Class 1 number of distinct alleles
  * 4-digit
	egen hla_class1_rownonmiss_l = rownonmiss(a_1-c_2), strok
	egen hla_class1_homo_l = rowtotal(hla_a_l_* hla_b_l_* hla_c_l_*) if hla_class1_rownonmiss_l==6
	* Note HLA_C 18 undifferentiated all have two C alleles, so can be included
	replace hla_class1_homo_l = hla_class1_homo_l+1 if c_2=="18"
	* Note HLA_A 74 undifferentiated has 2 A alleles, so can be included
	replace hla_class1_homo_l = hla_class1_homo_l+1 if a_2=="74"
	* Note 2 with HLA_C 17 undifferentiated are potentially homozygous - exclude 
	replace hla_class1_homo_l = hla_class1_homo_l+1 if c_1=="17"
	replace hla_class1_homo_l = hla_class1_homo_l+1 if c_2=="17"
	replace hla_class1_homo_l=. if c_1=="17" & c_2=="17"
	* Subtract 1 if allele is null (George said to send to missing - Perri agrees they have one less functional allele)
	replace hla_class1_homo_l = hla_class1_homo_l - 1 if hla_a_1_n==1
	replace hla_class1_homo_l = hla_class1_homo_l - 1 if hla_a_2_n==1
	
	* Recode to collapse categories
	recode hla_class1_homo_l 3/4=0 5=1 6=2, gen(hla_class1_homo_l_cat)
	label define hla_class1_homo_l_cat 0 "3/4" 1 "5" 2 "6"
	label val hla_class1_homo_l_cat hla_class1_homo_l_cat
	
	* Recode to make binary
	recode hla_class1_homo_l 3/5=0 6=1, gen(hla_class1_homo_l_b)
	label define hla_class1_homo_l_b 0 "3/5" 1 "6"
	label val hla_class1_homo_l_b hla_class1_homo_l_b
	
  * 2-digit
	egen hla_class1_rownonmiss_s = rownonmiss(hla_a_1_s hla_a_2_s hla_b_1_s hla_b_2_s hla_c_1_s hla_c_2_s), strok
	egen hla_class1_homo_s = rowtotal(hla_a_s_* hla_b_s_* hla_c_s_*) if hla_class1_rownonmiss_s==6
	
	* Missing if null allele - replace to be the same as homo_l variable - checked manually
	replace hla_class1_homo_s = hla_class1_homo_l if hla_class1_homo_s==.
	replace hla_class1_homo_s = 4 if pid=="MF2014" // manual count
	
	* Recode to collapse categories
	recode hla_class1_homo_s 3/4=0 5=1 6=2, gen(hla_class1_homo_s_cat)
	label define hla_class1_homo_s_cat 0 "3/4" 1 "5" 2 "6"
	label val hla_class1_homo_s_cat hla_class1_homo_s_cat
	
	* Recode to make binary
	recode hla_class1_homo_s 3/5=0 6=1, gen(hla_class1_homo_s_b)
	label define hla_class1_homo_s_b 0 "3/5" 1 "6"
	label val hla_class1_homo_s_b hla_class1_homo_s_b
	
 ** Class II number of distinct alleles
  * 4-digit
	egen hla_class2_rownonmiss_l = rownonmiss(dpa1_1-drb1_2), strok
	egen hla_class2_homo_l = rowtotal(hla_dpa1_l_* hla_dpb1_l_* hla_dqa1_l_* hla_dqb1_l_* hla_drb1_l_*) if hla_class2_rownonmiss_l==10
	
	* Recode to collapse categories
	recode hla_class2_homo_l 4/7=0 8=1 9=2 10=3, gen(hla_class2_homo_l_cat)
	label define hla_class2_homo_l_cat 0 "4/7" 1 "8" 2 "9" 3 "10"
	label val hla_class2_homo_l_cat hla_class2_homo_l_cat
	
	* Recode to make binary
	recode hla_class2_homo_l 4/9=0 10=1, gen(hla_class2_homo_l_b)
	label define hla_class2_homo_l_b 0 "4/9" 1 "10"
	label val hla_class2_homo_l_b hla_class2_homo_l_b
		
********************************************************************************	
********************************************************************************
	
*** SAVE DATA   
	drop *_length /* Drop extra vars */
	
	replace pid = substr(pid, 3, 6)
	rename pid id
	destring id, replace
	
 ** Order variables and save data	
	order _all, sequential
	order id, first
	
	save "`datapath'/`outputfile'", replace		

********************************************************************************
********************************************************************************

*** COHORT DATA

*** INDIVIDUAL - ONE RECORD PER INDIVIDUAL
	use "`datapath'/`inputfile_ind'", clear
	rename agecat agecat_enroll
	rename childadult child
	
 ** Male
	gen male=gender==1
	
 ** G6PD
	recode G6PD (0 1 3 = 0) (2 4 = 1) (9 = .), gen(g6pd_r)
	tab G6PD g6pd_r, m
	label var g6pd_r "G6PD: Male hemizygote/Female homozygote"
	
 ** HBS
	recode hbs (1 2= 1) (9 = .), gen(hbs_r)
	label var hbs_r "Hb AS/Hb SS"
	
 ** ALPHA THALASSEMIA
	recode alphathal (9 = .), gen(alphathal_r)
	label var alphathal_r "Alpha globin variant"
	label val alphathal_r alphathal

*** ALL STUDY VISITS - LONG FORM - MULTIPLE RECORDS PER INDIVIDUAL
	merge 1:m id using "`datapath'/`inputfile_clin'"
	drop _m
	
 ** DROP VISITS FOR TORORO (site 3) AFTER DEC 31, 2014
	drop if siteid==3 & date > date("20141231","YMD") /*CHECK DATE*/

*** INCIDENCE
 ** MALARIA INCIDENCE BY PERSON (ENTIRE STUDY)
  * PERSON-TIME FOR ENTIRE STUDY
	bysort id: egen firstvisit = min(date)
	bysort id: egen lastvisit = max(date)
	format %td firstvisit lastvisit
	gen pt_study = lastvisit-firstvisit
	browse id date firstvisit lastvisit pt_study
	
	gen pt_study_month = pt_study/(365/12)
	
  * EPISODES OF INCIDENT MALARIA
	bysort id: egen incidentmalaria_tot = total(incidentmalaria)
	
  * INCIDENCE RATE
	gen incidencerate = incidentmalaria_tot/pt_study
	
	browse id age date firstvisit lastvisit incidentmalaria* pt_study incidencerate
	
 ** MALARIA INCIDENCE BY PERSON BY YEAR
  * PERSON-TIME BY YEAR
   	* CREATE YEAR VARIABLE FOR EACH VISIT
	  gen year = year(date)
  
    * YEAR OF LAST VISIT
	  gen firstvisit_year = year(firstvisit)
	  
	* YEAR OF LAST VISIT
	  gen lastvisit_year = year(lastvisit)
	
	* PERSON-TIME PER YEAR
	  * Long form
	  gen pt_year = .
	  forvalues i = 2011/2016 {
		* for the first year a patient is in the study PT = year/12/31-first study visit
		replace pt_year = date("`i'1231","YMD")-firstvisit if firstvisit_year==`i' & year==`i'
		* for all in between years PT = 365 (except 2012 PT = 366: leap year)
		replace pt_year = 365 if firstvisit_year<`i' & lastvisit_year>`i' & year==`i'
		* for the last year a patient is in the study PT = last study visit - year/1/1
		replace pt_year = lastvisit-date("`i'0101","YMD") if lastvisit_year==`i' & year==`i'
		
		* add 1 to PT for 2012 due to leap year
		if `i' == 2012 { /* leap year */
			replace pt_year = pt_year + 1 if firstvisit_year<`i' & lastvisit_year>`i' & year==`i'
		}
	  }
	  
	  replace pt_year = . if firstvisit==lastvisit
	  replace pt_year = lastvisit-firstvisit if firstvisit_year==lastvisit_year
	  browse id date firstvisit lastvisit pt_*

	* EPISODES OF INCIDENT MALARIA BY YEAR
	  bysort id year: egen incidentmalaria_year = total(incidentmalaria)
	  browse id date incidentmalaria_year incidentmalaria
	
	* INCIDENCE RATE BY YEAR
	  gen incidencerate_year = incidentmalaria_year/pt_year
	  browse id age date firstvisit lastvisit year incidentmalaria_year pt_year incidencerate_year
	
	* CHECKS
	  /* keep id year pt_year incidencerate_year
		 duplicates drop
		 sum pt_year
		 sum incidencerate_year */

*** PARISITEMIA PREVALENCE (3 MONTH INTERVALS - BINARY Y/N)
	browse id date routinevisit age parasitedensity monthyear LAMP malariacat
	
 ** CREATE VISIT PERIOD
  * MIN DATE - 8/5/11
  * MAX DATE - 6/30/16
	gen month = month(date)

	gen q = 1 if month==1|month==2|month==3
		replace q = 2 if month==4|month==5|month==6
		replace q = 3 if month==7|month==8|month==9
		replace q = 4 if month==10|month==11|month==12
	
	bysort id year q: egen para_q_max = max(parasitedensity)
	bysort id year q: egen lamp_q = max(LAMP)
	
	recode para_q_max (0=0) (1/max=1), gen(para_q)
	
	gen para_lamp_q = 1 if (para_q==1|lamp_q==1)
		replace para_lamp_q = 0 if para_q==0 & (lamp_q==0|lamp_q==.)
		
	gen para_lamp_q_all = para_lamp_q /* LAMP including adults */
	
	replace para_lamp_q = . if child==0 /* Original LAMP variable for children only */
		/* LAMP not done on all adults */
		/* For now, there are 862 routine child visits for which parasite density = 0 and LAMP is missing,
		   7% of visits that should have had a LAMP;
		   include these for now - this is a limitation of the analysis */
	
	/* 	Send variable including adults to missing after dates for which LAMP on 
		adults stopped at each site because then adults only have microscopy data */
		forvalues i=1/3 {
			sum date if siteid==`i' & LAMP<. & child==0
			di %td r(max)
			replace para_lamp_q_all=. if siteid==`i' & date>r(max) & child==0
			
			sum date if siteid==`i' & para_lamp_q_all<. & child==0
			di %td r(max)
			
			sum date if siteid==`i' & para_lamp_q_all<. & child==1
			di %td r(max)
		}
		
 ** Person-time LAMP all (including children and adults)
	sort id date
	gen pt_study_lampall_lastvisit_raw = date if para_lamp_q_all<. 
	by id: egen pt_study_lampall_lastvisit = max(pt_study_lampall_lastvisit_raw)
	format %td pt_study_lampall_lastvisit*
	
	gen pt_study_lampall = pt_study_lampall_lastvisit-firstvisit
	gen pt_study_month_lampall = pt_study_lampall/(365/12)
	
	browse id child date firstvisit lastvisit pt_study pt_study_*lampall*
		
*** LOG-TRANSFORMED PARASITEMIA
	gen para_ln = log(parasitedensity)

*** AGE 
 ** AGECAT
	rename agecat agecat_visit
	
 ** MEAN AGE FOR YEAR
	bysort id: egen age_first = min(age)
	bysort id: egen age_last = max(age)

	gen age_year_start = .
	gen age_year_end = .
	
	forvalues i = 2011/2016 {
		/* Calculate age as the midpoint of the year
           Midpoint of entire year if not first or last year observed
           Midpoint of time from first visit to Dec 31 for first year
           Midpoint of time from Jan 1 to last visit for last year */
		* for all years except first, make start age age on Jan 1
		replace age_year_start = (date("`i'0101","YMD")-dob)/365.25 if firstvisit_year!=`i' & year==`i'
		* for first year, make start age first recorded age
		replace age_year_start = age_first if firstvisit_year==`i' & year==`i'
		* for all years except last, make end age age on Dec 31
		replace age_year_end = (date("`i'1231","YMD")-dob)/365.25 if lastvisit_year!=`i' & year==`i'
		* for last year, make end age last recorded age
		replace age_year_end = age_last if lastvisit_year==`i' & year==`i'
		}
	
	gen age_y = (age_year_end+age_year_start)/2
	
	gen agecat_y = 1 if age_y<5
		replace agecat_y = 2 if age_y>=5 & age_y<12
		replace agecat_y = 3 if age_y>=18 & age_y<.
	tab age_y if agecat_y==.
	label val agecat_y agecat
	
 ** MEAN AGE PER QUARTER
	gen age_q_start = .
	gen age_q_end = .
	
	bysort id: gen firstvisit_q_raw = q if _n==1
	bysort id: egen firstvisit_q = max(firstvisit_q_raw)
	bysort id: gen lastvisit_q_raw = q if _n==_N
	bysort id: egen lastvisit_q = max(lastvisit_q_raw)
	
	forvalues i = 2011/2016 {
			/* 	Calculate age as the midpoint of the quarter
				Midpoint of entire quarter if not first or last quarter observed
				Midpoint of time from first visit to end of quarter for first quarter
				Midpoint of time from start of quarter to last visit for last quarter */
			replace age_q_start = (date("`i'0101","YMD")-dob)/365.25 if year==`i' & q==1
			replace age_q_start = (date("`i'0401","YMD")-dob)/365.25 if year==`i' & q==2
			replace age_q_start = (date("`i'0701","YMD")-dob)/365.25 if year==`i' & q==3
			replace age_q_start = (date("`i'1001","YMD")-dob)/365.25 if year==`i' & q==4
			replace age_q_start = age_first if firstvisit_year==`i' & year==`i' & firstvisit_q==q
			
			replace age_q_end = (date("`i'0331","YMD")-dob)/365.25 if year==`i' & q==1
			replace age_q_end = (date("`i'0630","YMD")-dob)/365.25 if year==`i' & q==2
			replace age_q_end = (date("`i'0930","YMD")-dob)/365.25 if year==`i' & q==3
			replace age_q_end = (date("`i'1231","YMD")-dob)/365.25 if year==`i' & q==4
			replace age_q_end = age_last if lastvisit_year==`i' & year==`i' & lastvisit_q==q
	}
	
	gen age_q = (age_q_end+age_q_start)/2
	
	gen agecat_q = 1 if age_q<5
		replace agecat_q = 2 if age_q>=5 & age_q<12
		replace agecat_q = 3 if age_q>=18 & age_q<.
	tab age_q if agecat_q==.
	label val agecat_q agecat
	
	browse id age year q age_q*
	
*** MERGE IN KIR/HLA DATASET
	merge m:1 id using "`datapath'/`outputfile'"
	rename _m merge_clinical
	
 ** MISSING FROM KIR/HLA DATA		
	quietly tab id if merge_clinical==1
	display r(r)
	*486 in clinical data and NOT in KIR/HLA
	quietly tab id if merge_clinical==2
	display r(r)	
	*1 in KIR/HLA and NOT in clinical data

	keep if merge_clinical==3
	
	save "`datapath'/`outputfile'", replace

********************************************************************************

*** EIR: Household level
	use "`datapath'/`inputfile_EIR'", clear 

 ** Exclude time not in analysis
  * DROP VISITS FOR TORORO (site 3) AFTER DEC 31, 2014
	drop if siteid==3 & date > date("20141231","YMD") 

 ** Calculate P. falciparum Sporozoite rates by site
	bysort siteid: egen numberpositive_tot = total(numberpositive)
	bysort siteid: egen numbertested_tot = total(numbertested)
	gen pfsr = numberpositive_tot/numbertested_tot

 ** Calculate human biting rate: geometric average adding 0.5 to all values
	levelsof hhid, local(levels)
	
	gen hbr = .
	
	* Loop through all values of HHID /* Ameans doesn't work with egen */
	foreach l of local levels {
		ameans(totalanopheles) if hhid==`l', add(0.5)
		replace hbr = r(mean_g) if hhid==`l'
	}
	
 ** Calculate EIR: pfsr * hbr * 365 days
	gen eir = pfsr * hbr * 365
	
 ** Calculate log eir
	gen eir_ln = log(eir)
	
	keep hhid eir eir_ln
	duplicates drop
	
*** MERGE IN KIR/HLA DATASET
	merge 1:m hhid using "`datapath'/`outputfile'"
	rename _m merge_eir
	
	order id hhid siteid eir eir_ln

 ** MISSING FROM EIR DATA	
 	quietly tab id if merge_eir==3
	display r(r)
	*892 subjects are in EIR data, two IDs (1 HH) are missing EIR because it wasn't collected
	
	drop if hhid==201031206 /* Missing EIR */
	
	keep if merge_eir==3 /* EIR DATA COMPLETE */
	
	save "`datapath'/`outputfile'", replace
	
********************************************************************************

*** Language	
 ** Jinja
	import excel "`datapath'/`inputfile_langJIN'", sheet("Jinja cohort") firstrow case(lower) clear
	keep id language
	drop if id==1114 /* Inconsistent duplicate ID not in KIR data */
	
	tempfile languagedata    /* create a temporary file */
	save "`languagedata'"      /* save memory into the temporary file */
	
 ** Kanungu
	import excel "`datapath'/`inputfile_langKAN'", sheet("Kanungu") firstrow case(lower) clear
	keep id language
	
	append using "`languagedata'"
	save "`languagedata'", replace
	
 ** Kanungu
	import excel "`datapath'/`inputfile_langTOR'", sheet("Tororo cohort") firstrow case(lower) clear
	keep id language
	
	append using "`languagedata'"
	
	drop if id==3194 & language=="Muganda" /* Other HH members speak Jap */
			
 ** Merge to cohort data	
	merge 1:m id using "`datapath'/`outputfile'"
	keep if _m!=1 /* Drop if not in KIR data */
	save "`datapath'/`outputfile'", replace
	
 ** Create clean language variable
	gen language_r = lower(language)
	
  * Fill in using HH ID
	bysort hhid: egen languagemode = mode(language_r)
	replace language_r = languagemode if language_r==""
	drop languagemode
	
  * Clean language string
	replace language_r =  "luganda" if language_r == "ganda"
	replace language_r =  "ateso" if language_r == "ites0"
	replace language_r =  "ateso" if language_r == "iteso"
	replace language_r =  "ateso" if language_r == "itesot"
	replace language_r =  "dhopadhola " if language_r == "jaluo"
	replace language_r =  "dhopadhola " if language_r == "jap"
	replace language_r =  "rukiga" if language_r == "kiga"
	replace language_r =  "dhopadhola " if language_r == "ludama"
	replace language_r =  "rukiga" if language_r == "lukiga"
	replace language_r =  "runyankole" if language_r == "lunyankole"
	replace language_r =  "kinyarwanda" if language_r == "lunyarwanda"
	replace language_r =  "runyoro" if language_r == "lunyoro"
	replace language_r =  "samia" if language_r == "lusamya"
	replace language_r =  "ateso" if language_r == "luteso"
	replace language_r =  "lunyole" if language_r == "munyole"
	replace language_r =  "runyankole" if language_r == "nkore"
	replace language_r =  "kinyarwanda" if language_r == "rwandese"

 ** Ethnicity: 3 Category
	gen ethnicity = .
	replace ethnicity =  3 if language_r =="acholi"
	replace ethnicity =  2 if language_r =="ateso"
	replace ethnicity =  3 if language_r =="dhopadhola "
	replace ethnicity =  2 if language_r =="karamajong"
	replace ethnicity =  1 if language_r =="kikuyu"
	replace ethnicity =  1 if language_r =="kinyarwanda"
	replace ethnicity =  2 if language_r =="langi"
	replace ethnicity =  1 if language_r =="luganda"
	replace ethnicity =  3 if language_r =="lugbara"
	replace ethnicity =  1 if language_r =="lugisu"
	replace ethnicity =  1 if language_r =="lugwere"
	replace ethnicity =  1 if language_r =="lunyole"
	replace ethnicity =  1 if language_r =="lusoga"
	replace ethnicity =  1 if language_r =="lutoro"
	replace ethnicity =  1 if language_r =="rukiga"
	replace ethnicity =  1 if language_r =="runyankole"
	replace ethnicity =  1 if language_r =="runyoro"
	replace ethnicity =  1 if language_r =="samia"
	label define ethnicity 1 "Bantu" 2 "Nilo-Hamite" 3 "Luo"
	label val ethnicity ethnicity
	
 ** Bantu vs. non-Bantu
	recode ethnicity (2/3 = 0), gen(bantu)
	
 ** Site/Bantu
	gen site_bantu = siteid
		replace site_bantu = 0 if siteid==1 & bantu==1 
		replace site_bantu = . if siteid==1 & ethnicity==.
		label define site_bantu 0 "Jinja/Bantu" 1 "Jinja/Non-Bantu" 2 "Kanungu (99% Bantu)" 3 "Tororo (99% Non-Bantu)"
		label val site_bantu site_bantu
		
********************************************************************************
********************************************************************************		
		
*** Keep only necessary variables
	keep 	id date hhid siteid site_bantu bantu male eir_ln child g6pd_r hbs hbs_r alphathal_r ///
			pt_study_month pt_study_month_lampall incidencerate ///
			para_ln visittype ///
			para_q para_lamp_q para_lamp_q_all year q ///
			incidentmalaria_year pt_year ///
			parasitedensity temperature ///
			age age_q agecat_q age_y agecat_y ageyrs ///
			hla* febrile	
		
	save "`datapath'/`outputfile'", replace
	
********************************************************************************
********************************************************************************

* Probability of symptoms if infected dataset

*** Tag parasitemia to exclude from probability of symptoms if infected analysis	
 ** Exclude asymptomatic parasitemia where fever develops within 7 days (Keep the febrile parasitemia)
 ** Exclude febrile parasitemia where had febrile parasitemia in prior 7 days (Keep first febrile parasitemia only) 
	
	keep if parasitedensity > 0 & parasitedensity < . /* Keep visits that were parasitemic */
	sort id date

	* Tag aymptomatic parasitemia where fever develops within 7 days
	gen symp7days1 = 1 if (date[_n+1] - date) <=7 & parasitedensity > 0 & parasitedensity < . & febrile==0 & febrile[_n+1]==1 & id==id[_n+1]
	replace symp7days1 = 1 if (date[_n+2] - date) <=7 & parasitedensity > 0 & parasitedensity < . & febrile==0 & febrile[_n+2]==1 & id==id[_n+2]
	replace symp7days1 = 1 if (date[_n+3] - date) <=7 & parasitedensity > 0 & parasitedensity < . & febrile==0 & febrile[_n+3]==1 & id==id[_n+3]

	bysort id: egen symp7days1_max = max(symp7days1)
	browse id date symp7days* parasitedensity febrile visittype incidentmalaria if symp7days1_max==1
	
	* Tag febrile parasitemia where had prior febrile parasitemia within 7 prior days
	  * Tags the second instance of parasitemia to drop
	gen symp7days4 = 1 if date - date[_n-1] <=7 & parasitedensity[_n-1] > 0 & parasitedensity[_n-1] < . & febrile[_n-1]==1 & ///
												  parasitedensity > 0 & parasitedensity < . & febrile==1 & id==id[_n-1]
	browse id date symp7days* parasitedensity febrile visittype incidentmalaria if symp7days4==1

	* Tag parasitemia episodes to EXCLUDE from analysis
	gen exclude_pd = 1 if symp7days1==1|symp7days4==1
	
	* Drop parasitemia episodes to exclude and non-parasitemic visits
	drop if exclude_pd == 1|parasitedensity==0
	
	drop symp7days* exclude_pd date
	
	*** Save probability of infections	
	save "`datapath'/`outputfile_psi'", replace
	
********************************************************************************

*** Drop date
	use "`datapath'/`outputfile'", clear
	drop date
	save "`datapath'/`outputfile'", replace
	
