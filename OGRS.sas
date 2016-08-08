* OGRS for SPSS Version 1.2; 
* Copyright 2016;
* by Amanda Kay Montoya;
* Documentation available at www.afhayes.com;

* Permission is hereby granted, free of charge, to any person obtaining a copy of this software;
* and associated documentation files (the "Software"), to use the software in this form.  Distribution;
* after modification is prohibited, as is its use for any commercial purpose without authorization;  
* This software should not be posted or stored on any webpage, server, or directory accessible to;
* the public whether free or for a charge unless written permission has been granted by the copyright;
* holder.  The copyright holder requests that this software be distributed by directing users to;
* afhayes.com where the latest release of the software and related documentation is archived and;
* can be downloaded;

* THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, ;
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF ;
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT ;
* IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, ;
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT ;
* OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE ;
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE ;

* The above text should be included in any distribution of the software;

* Changes to this version include added error messages not present in Version 1.1;
	* CONF must be greater than or equal 50 and less than 100; 
	* Only one variable can be specified as the X argument;
    * Only one variable can be specified as the M argument;
    * Only one variable can be specified as the Y argument;
	* Missing data is deleted listwise;

%macro RedR (center = );
	mcent = m - &center;
	m2int = mcent#data[,2:numgroup];
	datam2 = J(Ninit,1,1)||mcent||m2int;
	IF (covtog = 1) then datam2 = datam2||cov;
	yestm2 = datam2*inv((datam2`)*datam2)*datam2`*y;
	ycorm2 = CORR(yestm2||y);
	ycorm2 = ycorm2[1,2];
	redr2 = ycorm2**2;
%mend;

%macro PROBE (min = , max = );
*transtog = 1;

*Space between probed points;
jump = (&max - &min)/iter;

dim = (transtog = 0)*(iter+1) + (transtog = 1)*(iter-1);

*Results columns: (1) M-value, (2) R^2 reduced, (3) Change in R^2;
* (4) F-stat test of change R^2, (5) p-value change in R^2, (6) Increasing/Decreasing Flag;
* (7) Transition above and below crit F, (8) Within convergence criteria;

tempres = J(dim, 8, -999);

IF (transtog = 0) then i2 = &min;
IF (transtog = 1) then i2 = &min + jump;

do i = 1 TO dim;
	%RedR(center = i2);
	tempres[i, 1:2] = i2||redr2;
	i2 = i2 + jump;
end; 

*Change in R^2 reduced to full;
tempres[,3] = fullr2 - tempres[,2];
*F statistic for change in R^2;
tempres[,4] = (dffull*tempres[,3])/(dfred*(1-fullr2));
*p-value for change in F;
tempres[,5] = 1 - CDF('F',tempres[,4],dfred,dffull);




%mend;

%macro OGRS (data =, vars =, x = , m = , y = , conf = 95, convcrit = .00000001, decimals = 10.4, iter = 0);
proc iml; 
runnotes = j(100,1,0);
runcount = 1;
criterr = 0;
use &data;
read all var{&vars} into allvars; 
allnames = {&vars};
read all var{&x} into xdat;
xname = {&x};
if ncol(xdat) > 1 then do;
	runnotes[runcount,1] = 2;
	runcount = runcount +1;
	criterr = 1;
end;
read all var{&m} into mdat;
mname = {&m};
if ncol(mdat) > 1 then do;
	runnotes[runcount,1] = 3;
	runcount = runcount +1;
	criterr = 1;
end;
read all var{&y} into ydat;
yname = {&y};
if ncol(ydat) > 1 then do;
	runnotes[runcount,1] = 4;
	runcount = runcount +1;
	criterr = 1;
end;

xx=(allvars=.);
xx=xx[,+];
j=1;do i=1 to nrow(allvars);if xx[i,1]=0 then;do;allvars[j,]=allvars[i,];j=j+1;end;end;
allvars=allvars[1:j-1,];
if sum(xx) > 0 then; do; runnotes[runcount,1] = 1; runcount = runcount +1; end;
missing = nrow(xx) - nrow(allvars);

if (criterr = 0) then;do;
  do i = 1 to ncol(allvars)-1;
  	do j = i+1 to ncol(allvars);
    	if (allnames[1,i]=allnames[1,j]) then;do;
      	runnotes[runcount,1]=8;
		runcount = runcount +1;
		criterr=1;
		end;
    end;
  end;
end;

if ((xname = mname)|(xname = yname)|(mname= yname)) then do; 
	runnotes[runcount,1] = 8; 
	runcount = runcount + 1;
	criterr = 1;
end;

zero=j(nrow(allvars),1,0);
do i = 1 to (ncol(allvars)-1);
  do j = i+1 to ncol(allvars);
  	diff=allvars[,i]-allvars[,j];
  	copy=(diff=zero);copy=copy[+];
  if (copy=nrow(allvars)) then;do;
    copyname=allnames[1,i]||allnames[1,j];
	runnotes[runcount,1]=9;
	criterr=1;
	runcount = runcount + 1;
  end;
  end;
end;

convcrit = &convcrit;
if (convcrit < .00000001) then; do; 
	convcrit = .00000001;
	runnotes[runcount,1] = 6;
	runcount = runcount+1;
end;
conf = &conf;
if ((conf < 50) | (conf > 99.9999)) then do;
	conf = 95;
	runnotes[runcount,1] = 5;
	runcount = runcount+1;
end;
alpha = 1-conf/100;

*Toggle for if there are covariates;
covtog = (ncol(allvars) - 3 > 0);

Ninit = nrow(allvars);

if covtog=1 then do;
	covcount = 1;
	cov = J(Ninit, ncol(allvars)-3, 999);
	covname = J(1, ncol(allvars)-3, "ERROR");
end;



do i = 1 TO ncol(allnames);
	if (allnames[,i] = xname) then do;
		x = allvars[,i];
		end;
	if (allnames[,i] = mname) then do;
		m = allvars[,i];
		end;
	if (allnames[,i] = yname) then do;
		y = allvars[,i];
		end;
	if  all(allnames[,i] ^= xname||mname||yname) then do;
		if covtog = 1 then do;
			cov[,covcount] = allvars[,i];
			covname[,covcount] = allnames[,i];
			covcount = covcount + 1;
		end;
	end;
end; 


designx = design(x);
numgroup = ncol(designx);
designx = designx[,1:(numgroup-1)];

if (criterr <> 1) then do;
	xmat = J(ncol(designx)+1, ncol(designx) + 1, -999);
	DO kloop = 1 TO ncol(designx);
		icount = 1;
		DO WHILE (xmat[kloop,1] = -999); 
			IF designx[icount,kloop] = 1 then xmat[kloop,1] = x[icount,1];
			icount = icount + 1;
		END;
	END;
	icount = 1;
	DO WHILE (xmat[ncol(designx)+1,1] = -999);
		IF all(designx[icount,] = 0) then do; 
			xmat[ncol(designx)+1,] = x[icount,1]||J(1,ncol(designX),0);
		end;
		icount = icount +1;
	END;
	xmat[1:(numgroup-1),2:numgroup] = I(numgroup-1);

	prodcol = designX#m;

	if covtog = 0 then data = J(Ninit,1,1)||designX||m||prodcol;
	if covtog = 1 then data = J(Ninit,1,1)||designX||m||prodcol||cov;

	*Calculate estimated outcome variable;
	yest = data*inv(data`*data)*data`*y;

	*Calculate R^2 for full model with all predictors;
	ycor = corr(y||yest);
	ycor = ycor[1,2];
	fullr2 = ycor**2;

	*NEED TO INSERT SOME CODE WHICH DETERMINES ITER;
	if (&iter ^=0) then;do; 
		iter=abs(round(&iter));
		if (&iter ^= iter) then do;
			runnotes[runcount, 1] = 7;
			runcount = runcount +1;
		end;
	end;
	if &iter = 0 then iter = 50+10*numgroup;

	*Degrees of freedom for full model. Participants minus predictors; 
	dffull = Ninit - ncol(data);
	*Degrees of freedom change to reduced model. Number of groups minus 1;
	dfred = numgroup - 1;
	Ffull = (fullr2*dffull)/((1-fullr2)*(ncol(data)-1));
	pfull = 1 - CDF('F',Ffull,(ncol(data)-1),dffull);
	*Critical F Statistic;
	critF = FINV(conf/100, dfred, dffull);

	*Regression results matrix;
	modres = J(ncol(data), 6, -999);
	modres[,1] = inv(data`*data)*data`*y;
	ssr = sum((y-yest)##2);
	msr = ssr/(Ninit - ncol(data));
	semat = msr*inv(data`*data); 
	modres[,2] = (vecdiag(semat))##(1/2);
	modres[,3] = modres[,1]/modres[,2];
	modres[,4] = 2*(1-CDF('t', abs(modres[,3]),dffull));
	tcrit = TINV(1-alpha/2, dffull);
	modres[,5] = modres[,1] - tcrit*modres[,2];
	modres[,6] = modres[,1] + tcrit*modres[,2];

	*Results for omnibus test of interaction;
	dataint = J(Ninit, 1,1)||designx||m;
	IF (covtog = 1) then dataint = dataint||cov;
	yestint = dataint*inv(dataint`*dataint)*dataint`*y;
	ycorint = CORR(yestint||y);
	ycorint = ycorint[1,2];
	r2int = ycorint##2;
	rchint = fullr2 - r2int;
	Fint = (dffull*rchint)/(dfred*(1-fullr2));
	pint = 1 - CDF('F', Fint, dfred, dffull);
	intres = rchint||Fint||dfred||dffull||pint;



	*Transition toggle, will turn to 1 when a certain point is deemed a "transition";
	transtog = 0;
	*Calculate minimum and maximum values of the moderator;
	minM = min(m);
	maxM = max(m);

	%PROBE (min = minM, max = maxM);

	results = tempres;
	OGres = tempres;

	*Decreasing/Increasing flat and transition flag for last row;
	results[nrow(results),6:7] = {0 0};

	i3 = 1;

	DO WHILE (i3 <= nrow(results));
		IF(i3 < nrow(results)) then; do;
			*Tag increasing and decreasing;
			results[i3,6] = (results[i3,4] < results[i3+1,4])-(results[i3,4] > results[i3+1,4]);
			*Tag transitions above and below critical F;
			*Negative transition point: This point above F and row below is below F;
			*Positive transition point: This point below F and row before is above F;
			results[i3,7] = -1*((results[i3,4] > critF) & (results[i3+1,4] < critF)) + ((results[i3,4] < critF) & (results[i3+1,4] > critF));
		END;
		*Tag JN points (close enough by convergence criteria);
		results[i3,8] = (abs(results[i3,4] - critF) < convcrit);
		IF (i3 = nrow(results)) then do;
			transcnv = 0;
		end;
		IF (i3 = 1) then do;
			transcnv = ((results[i3,7] = 1)&((results[i3,8]=1)|(abs(results[i3+1,8] - critF) < convcrit)));
		end;	
		IF ((i3 ^= nrow(results))&(i3 ^= 1))then do;
			*Point is positive transition point (increasing) and this point or the point in the below row converges;
			trnscnv1 = ((results[i3,7] = 1) & ((results[i3,8] = 1) | (abs(results[i3+1,8] - critF) < convcrit)));
			*Point is negative transition point and this point or the point in the below row converges;
			trnscnv2 = ((results[i3,7] = -1) & ((results[i3,8] = 1) | (abs(results[i3+1,8] - critF) < convcrit)));
			transcnv = ((trnscnv1 = 1) | (trnscnv2 = 1));
		end;
		IF ((abs(results[i3,7]) = 1) & (transcnv = 0)) then do;
			trnsindx = i3;
			transtog = 1;
			minmtran = min(results[i3+1,1]||results[i3,1]);
			maxmtran = max(results[i3+1,1]||results[i3,1]);
			%PROBE (min = minmtran, max = maxmtran);
			results = results//tempres;
			call sort(results,1);
		end;
		IF ((abs(results[i3,7]) = 0) | (transcnv = 1)) then i3 = i3+1;	
	END;

	*Number of JN solutions. 1) Check if last row is soln, 2)Check if first row is solution but not transition, 3) Count transitions;
	numJN = (results[nrow(results),8]=1)+((results[1,8] = 1) & (results[1,7] ^= 1)) + sum(abs(results[,7]));

	IF (numJN > 0)then do;
		*Make empty matrix to fill with solutions;
		JNSoln = J(numJN, 1, -999);
		JNIndx = J(numJN, 1, -999);
		*Soln counter;
		slncnt = 1;
		IF (results[nrow(results),8] = 1) then do;
			JNSoln[1,1] = results[nrow(results),1];
			JNIndx[1,1] = nrow(results);
			slncnt = slncnt + 1;
		end;
		DO i1 = 1 to nrow(results);
			IF (abs(results[i1,7]) = 1) then do;
				abvblw = (results[i1,1]||abs(results[i1,4] - critF))//(results[i1+1,1]||abs(results[i1+1,4] - critF));
				unsort = abvblw;
				call sort(abvblw,2);
				JNSoln[slncnt,1] = abvblw[1,1];
				indxtog = all(abvblw = unsort);
				IF (indxtog = 1) then JNIndx[slncnt,1] = i1;
				IF (indxtog = 0) then JNIndx[slncnt,1] = i1+1;
				slncnt = slncnt + 1;
			end;
		end;
	end;
end;

print "************************ OGRS Procedure for SAS Version 1.2 *************************";
print "Written by Amanda K. Montoya";
print "Documentation available at www.afhayes.com";
print "********************************************************************************************";

if (criterr <> 1) then; do;
	varrlabs = {"X = " "M = " "Y = "}; 
	print (xname//mname//yname) [label = "Variables:" rowname = varrlabs];
	IF (covtog = 1) then do;
		print covname [label = "Statistical Controls:"];
	END;
	dummylab = {"D1" "D2" "D3" "D4" "D5" "D6" "D7" "D8" "D9"};
	xmatlab = xname||dummylab[1,1:(numgroup-1)];
	print xmat [label = "Dummy Variable Coding Scheme:" colname = xmatlab];
	print Ninit [label = "Sample Size:"];
	print "********************************************************************************************";
	print yname [label = "Outcome:"];
	modsum = sqrt(fullr2)||fullr2||Ffull||(ncol(data)-1)||dffull||pfull;
	print modsum [label = "Model Summary" colname = {"R" "R-Sq" "F" "df1" "df2" "p"} format = &decimals];
	intlab = {"Int1" "Int2" "Int3" "Int4" "Int5" "Int6" "Int7" "Int8" "Int9"};
	modlabs = "Constant"||dummylab[1,1:(numgroup-1)]||mname||intlab[1,1:(numgroup-1)];
	IF (covtog = 1) then modlabs = modlabs||covname;
	print modres [label = "Model" rowname = modlabs colname = {"coeff" "SE" "t" "p" "LLCI" "ULCI"} format = &decimals];

	intmat = J((numgroup-1),5,"ERROR");
	intmat[,1] = (intlab[1,1:(numgroup-1)])`;
	intmat[,2] = J((numgroup-1),1, "=");
	intmat[,3] = (dummylab[1,1:(numgroup-1)])`;
	intmat[,4] = J((numgroup-1),1,"X");
	intmat[,5] = J((numgroup-1), 1, mname);
	print intmat [label = "Interactions:"];

	print intres [label = "R-Square increase due to interaction(s):" colname = {"R2-chng" "F" "df1" "df2" "p"} format = &decimals];

	print "************************* JOHNSON-NEYMAN TECHNIQUE *************************";

	IF (iter > 10) then do;
	last = nrow(OGres);
	rjump = ceil(last/20);
	rowsel = 1;
	rcount = 1+rjump;
	DO WHILE (rcount <= last);
		rowsel = rowsel||rcount;
		rcount = rcount +rjump;
	END;
	IF (rcount-rjump ^= last) then rowsel = rowsel||last;
	END;
	JNtabnam = mname||"R2-chng"||"F"||"p";

	IF (numJN > 0) then do;
		print JNSoln [label = "Moderator value(s) defining Johnson-Neyman boundaries of significance" format = &decimals];
		IF (iter > 10) then do;
			JNouttab = OGres[rowsel,]//results[JNIndx,];
		END;
		IF (iter <= 10) then do;
			JNouttab = OGres//results[JNIndx,];
		END;
		call sort(JNouttab,1);
		JNouttab = JNouttab[,{1 3 4 5}];
	END;
	IF (numJN = 0) then do;
		print "No Johnson-Neyman bounds found within the range of observed data";
		IF (iter > 10) then JNouttab = OGres[rowsel,{1 3 4 5}];
		IF (iter <= 10) then JNouttab = OGres[,{1 3 4 5}];
	END;
	print JNouttab [label = "Conditional effect of X on Y at values of the moderator (M)" colname = JNtabnam format = &decimals];

end;

print "*************************** ANALYSIS NOTES AND WARNINGS ***************************";
print "Check SAS log for errors.  Do not interpret output if errors are found.";
print (nrow(results)) [label = "Number of points probed in Johnson-Neyman algorithm"];
do i = 1 to nrow(runnotes);
  if (runnotes[i,1]=1) then;do;
    print missing [label="NOTE: Some cases were deleted due to missing data.  The number of cases was:"];
  end;
  if (runnotes[i,1]=2) then;do;
    print "ERROR: Only one X variable can be specified in the X= list.";
  end;
  if (runnotes[i,1]=3) then;do;
    print "ERROR: Only one M variable can be specified in the M= list.";
  end;
  if (runnotes[i,1]=4) then;do;
    print "ERROR: Only one Y variable can be specified in the Y= list.";
  end;
  if (runnotes[i,1]=5) then;do;
    print "NOTE: The confidence specified was not between 50 and 99.9999.";
	print "Level of confidence was adjusted to 95%";
  end;
  if (runnotes[i,1]=6) then;do;
    print "NOTE: Convergence criteria specified was not greater than or equal to .00000001.";
	print "Level of convergence criteria was adjusted to .00000001.";
  end;
  if (runnotes[i,1]=7) then;do;
    print "NOTE: An invalid number of initial iterations was specified.";
	print iter [label = "The number of initial iterations used was:"];
  end;
  if (runnotes[i,1]=8) then;do;
    print "ERROR: All variables must be unique. No variables can be the same among X, M, and Y";
  end;
  if (runnotes[i,1]=9) then;do;
    print copyname [label = "ERROR: Two of the specified variables are copies.  The variable names are:"];
  end;
end;



print "*****************************************************************************";

quit;
%mend;
