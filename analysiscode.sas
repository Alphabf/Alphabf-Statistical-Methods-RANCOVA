
proc datasets library=work memtype=data kill; quit;

/*********************/
data mybase;
	set sasdata.ads;
	where time=0;
run;
proc sort; by testcd phase  ;  run;
/*********************/
proc mixed data=mybase;
  		by testcd phase;
		class dose;
		model CHANGE = DOSE ;
        LSMEANS DOSE;
          ods output
          LSMeans = lsm_base
        ;
run;


data baselsm;
	set lsm_base;
    SUBPHASE=1;
run;


/**********************************/
data mydata;
	set sasdata.ads;
	where time^=0;
run;

proc sort; by testcd phase  ;  run;

%macro RANCOVA_PARALLEL(type=, random= );

proc mixed data=mydata;
  		by testcd phase ; 
   		class subjid dose subphase; 

        model change = base dose subphase dose*subphase; 
      	&random
		repeated SUBPHASE / type=&type sub=subjid;
		lsmeans DOSE*SUBPHASE ;
		%let type = %sysfunc(compress(%superQ(type),%str(%(%))));

        ods output
              Tests3 =&type.tests (where=(Effect="*SUBPHASE"))
              LSMeans = &type.lsm
			  FitStatistics =&type.Fit
        ;
run;

data &type.Fit (drop=Descr);
	set &type.Fit;
	length type $15;
	/*    where=(Descr="AIC (Smaller is Better)");*/
	where Descr="AICC (Smaller is Better)";
	type="&type.   ";
run;


data &type.tests ;
	set &type.tests;
	length type $15;
	type="&type.   ";
run;

data &type.lsm ;
	set &type.lsm;
	length type $15;
	type="&type.   ";
run;

%mend RANCOVA_PARALLEL;


/***************************************************************/

%RANCOVA_PARALLEL(type=CS)
%RANCOVA_PARALLEL(type=CSH)

%RANCOVA_PARALLEL(type=AR(1) , random=%str(random subjid;))
%RANCOVA_PARALLEL(type=ARH(1), random=%str(random subjid;) )
/****************** Checking the fit to select the best model****/
data fitdata;
	set Csfit Cshfit Ar1fit  Arh1fit;
run;

proc sort data=fitdata; by testcd phase value; run;

data bestfit;
	set fitdata;
    by testcd phase value;
	if first.phase;
run;
proc sort; by testcd phase type; run;

/*************Stack the data and select the best************/
data lsm; set Cslsm Cshlsm Ar1lsm Arh1lsm ; run;
proc sort; by testcd phase  type; run;

data lsm_s;
	merge lsm (in=a) bestfit (in=b);
	by testcd phase  type;
	if  a and b;
run;

data A_lsm_s;
	set lsm_s baselsm;
run;


proc sort data=A_lsm_s; by testcd   phase SUBPHASE; run;

data report_lsm_ (rename=(dose=trtn SUBPHASE=atptn));
	set A_lsm_s;
	length rownum 8 rowlbl $40 cell $20; 
	retain  rownum;
	by testcd   phase SUBPHASE;
	rownum =4;
		rowlbl = 'LSM';
		if abs(Estimate) < 10 then cell   = strip(put(round(Estimate,&round4), &STATFMT4D));
		else if 10 =< abs(Estimate) < 100 then cell   = strip(put(round(Estimate,&round3), &STATFMT3D));
		else if 100 =< abs(Estimate) < 1000 then cell   = strip(put(round(Estimate,&round2), &STATFMT2D));
		else if abs(Estimate) =< 1000 then cell   = strip(put(round(Estimate,&round1), &STATFMT1D));
	output;
	rownum +1;
		rowlbl = 'LSM s.e.';
		if abs(Estimate) < 10 then cell   = strip(put(round(StdErr,&round5), &STATFMT5D));
		else if 10 =< abs(Estimate) < 100 then cell   = strip(put(round(StdErr,&round4), &STATFMT4D));
		else if 100 =< abs(Estimate) < 1000 then cell   = strip(put(round(StdErr,&round3), &STATFMT3D));
		else if abs(Estimate) =< 1000 then cell   = strip(put(round(StdErr,&round2), &STATFMT2D));
	output;
	keep testcd PHASE DOSE SUBPHASE rownum rowlbl cell;
run;

data report_lsm;
	set report_lsm_;
	if missing(atptn) then atptn=999;
run;

/*******************************/
/*******************************/

proc sort data=A_lsm_s (rename=(dose=trtn SUBPHASE=atptn)) out=A_lsm_ss ; by testcd phase  atptn trtn ; run;

data report_lsmPv ;
	set A_lsm_ss;
	length rownum 8 rowlbl $40 cell $20; 
	retain  rownum;
	
	by testcd phase  atptn trtn ;
	rownum =6;
		rowlbl = 'LSM = 0 p-value     ';
		cell   = strip(strip(put(Probt,psignif.)));
	output;
	keep testcd phase  rownum rowlbl cell atptn trtn;
run;


/********************************************************/
*  Create the statistics for Section 1 of the report;
proc summary print data=sasdata.ads missing stackods 
             n mean std;
  class testcd  dose phase SUBPHASE;
  var CHANGE;
  ods output Summary=bySUBPHASE;
run; quit;

proc summary print data=sasdata.ads (where=(PHASE^=0)) missing stackods   
             n mean std;
  class subjid testcd  dose phase;
  var CHANGE;
  ods output Summary=tempsum;
run; quit;


proc summary print data=tempsum missing stackods    
             n mean std;
  class testcd  dose phase;
  var Mean;
  ods output Summary=overall;
run; quit;

data summary;
	set bySUBPHASE overall;
	if missing(SUBPHASE) then SUBPHASE=999;
run;

proc sort data=summary; by testcd  dose phase SUBPHASE; run;  

*;
*  Create a separate row from specific column values of the 
*  Summary data - Section 1.
*;

data report_summary_ (rename=(SUBPHASE=atptn dose=trtn));
	set work.summary;
	length rownum 8 rowlbl $40 cell $20; 
	retain  rownum;
	section = 1;
	by testcd  dose phase SUBPHASE;
	rownum =1;
		rowlbl = 'Mean Change from Vehicle';
		if abs(mean) < 10 then cell   = strip(put(round(mean,&round4), &STATFMT4D));
		else if 10 =< abs(mean) < 100 then cell   = strip(put(round(mean,&round3), &STATFMT3D));
		else if 100 =< abs(mean) < 1000 then cell   = strip(put(round(mean,&round2), &STATFMT2D));
		else if abs(mean) =< 1000 then cell   = strip(put(round(mean,&round1), &STATFMT1D));
	output;
	rownum +1;
		rowlbl = 'SD';
		if abs(mean) < 10 then cell   = strip(put(round(StdDev,&round5), &STATFMT5D));
		else if 10 =< abs(mean) < 100 then cell   = strip(put(round(StdDev,&round4), &STATFMT4D));
		else if 100 =< abs(mean) < 1000 then cell   = strip(put(round(StdDev,&round3), &STATFMT3D));
		else if abs(mean) =< 1000 then cell   = strip(put(round(StdDev,&round2), &STATFMT2D));
	output;

	rownum +1;
	rowlbl = 'N';
	cell   = strip(put(N, &STATFMT0D));
	output;

	keep testcd phase  rownum rowlbl cell SUBPHASE dose; ;
run;
data report_summary;
	set report_summary_;
	if atptn=999 and rownum=2 then cell="NA";
run;


/*****************Stack for report *********************/

data tempdata1 (where=(phase^=0));
	set Report_lsmPv Report_lsm report_summary;
	if missing(atptn) then atptn=999;
	paramn=input(testcd,param.);
run;

proc sort ; by testcd paramn phase  trtn rownum rowlbl atptn;run;


*  Transpose the data to get the layout needed for the report for post dose ;

proc transpose data=tempdata1 prefix=time_
               out=report_post (drop=_name_);
  by testcd paramn phase  trtn rownum rowlbl;
  var cell;
  id atptn;
run; 


/***********base****************/
data tempdata2 (where=(phase=0));
	set Report_lsmPv Report_lsm report_summary;
	paramn=input(testcd,param.);
	atptn=0;
run;

proc sort ; by testcd paramn phase  trtn rownum rowlbl atptn;run;


*  Transpose the data to get the layout needed for the report for pre dose;

proc transpose data=tempdata2 prefix=time_
               out=report_base_ (drop=_name_);
  by testcd paramn phase  trtn rownum rowlbl;
  var cell;
  id atptn;
run; 

data report_base;
	set report_base_;
	phase=1;
run;

/******all the table ***********************/

data sasdata.reportfinal;
	merge report_base report_post;
	by testcd paramn phase  trtn rownum rowlbl;
run;
 







