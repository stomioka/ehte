***************************************************************************************************************;
* This program is written by Jinglin Zhong for others outside SMPA to run the eHTE analysis                   *;
* This program creates two plots: culmulative1.png and ITE1.png, need to put your dir in two gpath statement  *;                                  
* Date: 2/7/2024                                                                                              *;   
* modified 2/12/2024                                                                                          *;
* Need to run eHTE_macro.sas before runing this program                                                       *;   
***************************************************************************************************************;

*****replace yourdir by the real path to the place you saved eHTE_macro.sas;
%include "yourdir\eHTE_macro.sas";
libname outdata "yourdir\";
options mprint mlogic;


******The input dataset for all the calculation into data one***;
******One patient per row only keep the LOCF at the endpoint for example Week 6;
******for each patient, only need two variables: chg(you endpoint which maybe change from baseline at week 6)
                                                 trt01pn(=1 for placebo, =2 for trt dose 1, =3 for trt dose 2);


data one;
  set outdata.PSIL201_chg;
  keep trt01pn chg;
run;

ods _all_ close;
ods html;
ods listing gpath="V:\380-000\bipolar\380-301\csr\stat\pgm\exploratory\JZ";
*****************plot the culmulative distribution by trt01pn;
*****************you can adjust the size of the plot in the statement below;
*******replace the gpath to your own directory to save the plots;
ods graphics /reset reset=index noborder width=5in height=6in;


ods graphics /imagename="Histgram" outputfmt=png;
%plothist(one);

ods graphics /imagename="Culmulative" outputfmt=png;
%plotcul(one);


*****************calculate the eHTE and p-value and plot the ITE;
*****************you can adjust the size of the plot in the statement below;
ods graphics /reset reset=index noborder width=2in height=6in   ; 
ods graphics /imagename="ITE" outputfmt=png;
 %eHTE_p(one);
