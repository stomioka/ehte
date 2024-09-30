*****************************************************************************************;
** This program contains several macro for the eHTE estimator and pvalue calculation*****;
** Written by: Jinglin Zhong                                                           **;
** Date: Nov 14, 2023																   **;	
** there are two different ways: 1. use the 48 pct from 3 to 97 by 2 for every study  **;
**                               2.use the real number of patients in the active arms  **;
** The research team decided to use method 1. Both ways are programmed here            **;
** The only change to switch to method 2 is go to the end of this program              **;
**  delete the star at the beginning of Lines 275-277 and add a star to Lines 280-282  **;
*****************************************************************************************;

********This macro can be used for both real dataset and simulated dataset;
********For real data set, leave nperms blank;
%macro sigma(indatap, indatat,nperms, outdata, outdata_ind);
proc univariate data=&indatap. noprint;
  by &nperms. trt01pn;
  var chg;
  output out=pct_pcb pctlpts=0 to 100 pctlpre=p;
run;
proc transpose data=pct_pcb out=pct_pcb_t(keep=&nperms. _name_ col1);
  by &nperms. trt01pn;
  var p0-p100;
run;
data pct_pcb_t1(rename=(col1=pcb_chg));
  set pct_pcb_t;
  percentile=1*compress(_name_,"p");
run;
*****calculate the percentile for each of the observations in each treatment arm;
*****use group=101 to make the number of pct matches the number of pct from proc univariate on pcb;
proc rank data=&indatat. out=pct_trt group=101;
  by &nperms. trt01pn;
  var chg;
  ranks percentile;
run;

proc sort data=pct_trt;
  by &nperms. percentile trt01pn;
run;

********using matched placebo percentile response calculate the treatment effect for each observation in treatment arms;
data matched;
  merge pct_trt(in=aa) pct_pcb_t1;
  by &nperms. percentile;
  if aa;
  ITE=chg-pcb_chg;
run;
proc sort data=matched;
  by &nperms. trt01pn;
run;

proc means data=matched noprint;
  by  &nperms. trt01pn;
  where percentile<=97 and percentile>=3;
  var ITE;
  output out=&outdata. stddev=sigma;
run;
data &outdata_ind.;
  set matched;
  where percentile<=97 and percentile>=3;
run;

/*proc datasets lib=work;
  delete matched pct_pcb pct_trt pct_pcb_t pct_pcb_t1;
run;
quit;*/
%mend;

%macro sigma_48(indatap, indatat,nperms, outdata, outdata_ind);
***calculated all the percentiles from 3rd to 97th based on placebo response and treatment response; 
proc univariate data=&indatap. noprint;
  by &nperms. trt01pn;
  var chg;
  output out=pct_PCB pctlpts=3 to 97 by 2 pctlpre=p;
run;
proc univariate data=&indatat. noprint;
  by &nperms. trt01pn;
  var chg;
  output out=pct_trt pctlpts=3 to 97 by 2 pctlpre=p;
run;

proc transpose data=pct_pcb out=pct_pcb_t(keep=&nperms. trt01pn _name_ col1);
  by &nperms. trt01pn;
run;
proc transpose data=pct_trt out=pct_trt_t(keep=&nperms. trt01pn _name_ col1);
  by &nperms. trt01pn;
run;
proc sort data=pct_pcb_t;
  by &nperms. _name_;
run;
proc sort data=pct_trt_t;
  by &nperms. _name_;
run;

data pct_all;
  merge pct_trt_t(in=aa) pct_pcb_t(drop=trt01pn rename=(col1=col_pcb));
  by &nperms. _name_;
  if aa;
  ITE=col1-col_pcb;
run;
proc sort data=pct_all;
  by trt01pn &nperms.;
run; 

***cacluate the stddev based on 48 percentiles regardless how many observations in the treatment arm;
proc means data=pct_all noprint;
  by trt01pn &nperms.;
  var ITE;
  output out=&outdata. stddev=sigma;
run;
data &outdata_ind.;
  set pct_all;
run;
/*proc datasets lib=work;
  delete pct_all pct_trt_t  pct_pcb pct_trt pct_pcb_t;
run;
quit;
*/
%mend;

%macro pvalue(method=);
proc sort data=SimeHTE&method.;
  by trt01pn;
run;
data count;
  merge eHTE&method.(keep=trt01pn sigma) simeHTE&method.(rename=(sigma=sim_sigma));
  by trt01pn;
 if sim_sigma>=sigma then ind=1;
  else ind=0;
run;
proc means data=count noprint;
  by trt01pn;
  var ind;
  output out=pvalue(keep=trt01pn  pvalue) mean=pvalue;
run;
  
data all;
  merge eHTE&method. pvalue;
  by trt01pn;
  eHTE=sigma/&s_1.;
run;

proc print data=all; 
  var trt01pn sigma eHTE pvalue;
  title "eHTE and p-value";
run;
title;
data ITE&method.p;
  set ITE&method.;
  percentile=Compress(_name_, "pP")*1;
run;
proc sort data=ITE&method.p;
  by trt01pn percentile;
run;
data ITE&method._PCB;
  trt01pn=1;
  do percentile=3 to 97 by 2;
  ITE=0;
  output;
  end;
run;
data ite&method.plot;
  set ITE&method._PCB ITE&method.p;
run;

ods html;
proc sgplot data=ITE&method.plot;
 title;
 styleattrs datacontrastcolors=(blue red orange);
  scatter x=ITE y=percentile/group=trt01pn markerattrs=(symbol=circle ) name="Treatment" legendlabel="Treatment" ;
  refline 0 /axis=x;
  keylegend "Treatment";
run;

%mend;

%macro plothist(indata);
proc sort data=&indata.;
  by trt01pn;
run;
proc sgplot data=&indata.;
  title "Overlapping Histograms of Sx Score";
styleattrs datacolors=(blue  red orange) datacontrastcolors=(blue  red orange);
  histogram chg/group=trt01pn transparency=0.5;
  density chg/group=trt01pn type=kernel lineattrs=(thickness=3);
  xaxis label="Sx Score";
  yaxis label="Frequency";
run;title;
%mend;


%macro plotcul(indata);
proc rank data=&indata. out=ranksplots groups=100;
  var chg;
  ranks percentile;
  by trt01pn;
run;
proc sort data=ranksplots;
  by trt01pn percentile;
run;
proc sgplot data=ranksplots;
  title "Culmulative Scores Dists";
styleattrs datacontrastcolors=(blue  red orange);
  series y=percentile x=chg/group=trt01pn;
  xaxis label="Sx Scores";
  yaxis label="culmulative percentile";
run; title;
%mend;


%macro eHTE_p(indata);
proc sort data=&indata.;
  by trt01pn;
run;
*****getting summary stat from the original dataset;
proc means data=&indata. noprint;
  by trt01pn;
  var chg;
  output out=sumstat(drop=_type_ rename=(_freq_=nobs)) mean=m stddev=s;
run;
****creat a macro variable Ntrt to store the number of treatment;
proc sql noprint;
  select count(distinct trt01pn) into :Ntrt
  from sumstat;
quit;
*****create macro varialbes to store the num of patients, mean and stddev for each trt01pn value;
data _null_;
  set sumstat;
  call symputx("m_"||trim(left(trt01pn)), m);
  call symputx("s_"||trim(left(trt01pn)), s);
  call symputx("nobs_"||trim(left(trt01pn)), nobs);
run;
******for eHTE calculation;
data placebo treatment;
  set &indata.;
  if trt01pn=1 then output placebo;
  else output treatment;
run;

********simulate data;

****create a separate data set for simulated placebo data because the procedures will be different from other arms;
****Simulate the placebo arm based on placebo mean and standard deviation;

data simu1;
  call streaminit(123);
  do nperms=1 to 10000;
  trt01pn=1;
  do pt=1 to &nobs_1.;
    chg=rand("Normal",&m_1.,&s_1. );
	output;
  end;
  end;
run;


****simulate data for other treatment arms with their observed mean but placebo standard deviation;
data simutrt;
  call streaminit(456);
  do nperms=1 to 10000;
    %do i=2 %to &Ntrt;
      do trt01pn=&i.;
        do pt=1 to &&nobs_&i.;
           m=&&m_&i.;
	       chg=rand("Normal",m,&s_1.);
	       output;
        end;
      end;
    %end;
  end;
run;



****calculate sigma for the real dataset and also the simulated dataset;
*%sigma(placebo, treatment, ,eHTE, ITE);
*%sigma(simu1, simutrt, nperms,SimeHTE,simITE);
*%pvalue(method=);

****calculate sigma based on 48 percentiles for the real dataset and also the simulated dataset;
%sigma_48(placebo, treatment, ,eHTE_48, ITE_48);
%sigma_48(simu1, simutrt, nperms,SimeHTE_48, simITE_48);
%pvalue(method=_48);

%mend;
  
