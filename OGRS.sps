/*OGRS for SPSS Version 1.2*/. 
/* Copyright 2016 */.
/* by Amanda Kay Montoya */.
/* Documentation available akmontoya.com*/.

/* Permission is hereby granted, free of charge, to any person obtaining a copy of this software*/.
/* and associated documentation files (the "Software"), to use the software in this form.  Distribution*/.
/* after modification is prohibited, as is its use for any commercial purpose without authorization*/.
/* This software should not be posted or stored on any webpage, server, or directory accessible to*/.
/* the public whether free or for a charge unless written permission has been granted by the copyright*/.
/* holder.  The copyright holder requests that this software be distributed by directing users to*/.
/* afhayes.com where the latest release of the software and related documentation is archived and*/.
/* can be downloaded*/.

/* THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,*/.
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF*/.
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT*/.
/* IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,*/.
/* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT*/.
/* OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE*/.
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE*/.

/* The above text should be included in any distribution of the software.*/.

/* Changes to this version include added error messages not present in Version 1.1. */.
         /* CONF must be greater than or equal to 50 and less than 100.*/.
         /* Only one variable can be specified as the X argument.*/.
         /* Only one variable can be specified as the M argument.*/.
         /* Only one variable can be specified as the Y argument.*/.
         /* Missing data is deleted listwise. */.
         /* No variables can be entered multiple times in the VARS list. */
         /* No variables in the vars list may be assigned to multiple roles (X, M, or Y). */

preserve. 
set printback=off.

define CDFINVT (p = !charend('/') /df = !charend('/')). 
compute p0=-.322232431088.
compute p1 = -1.
compute p2 = -.342242088547.
compute p3 = -.0204231210245.
compute p4 = -.0000453642210148.
compute q0 = .0993484626060.
compute q1 = .588581570495.
compute q2 = .531103462366.
compute q3 = .103537752850.
compute q4 = .0038560700634.
compute ppv = !p.
  do if (!p > .5).
    compute ppv = 1-!p.
  end if.
  compute y5=sqrt(-2*ln(ppv)).
  compute xp=y5+((((y5*p4+p3)*y5+p2)*y5+p1)*y5+p0)/((((y5*q4+q3)*y5+q2)*y5+q1)*y5+q0).
  do if (!p <= .5).
    compute xp = -xp.
  end if.
compute toutput = sqrt(!df*(exp((!df-(5/6))*(xp**2)/(!df-(2/3)+.1/!df)**2)-1)). 
!enddefine. 

define CORR (var1 = !charend('/') /var2 = !charend('/')). 
   COMPUTE var1 = !var1. 
   COMPUTE var2 = !var2. 
   COMPUTE MeanV1 = csum(var1)/nrow(var1). 
   COMPUTE MeanV2 = csum(var2)/nrow(var2). 
   COMPUTE Var1Cent = var1 - MeanV1. 
   COMPUTE Var2Cent = var2 - MeanV2. 
   COMPUTE crosprod = csum(Var1Cent &* Var2Cent). 
   COMPUTE Var1SS = csum(Var1Cent &* Var1Cent). 
   COMPUTE Var2SS = csum(Var2Cent &* Var2Cent).
   COMPUTE rPears = crosprod / (sqrt(var1SS)*sqrt(var2SS)). 
!enddefine. 

define RedR (center = !charend('/')). 
   COMPUTE mcent = m - !center. 
   COMPUTE m2int = MAKE(N, numgroup-1, -999).
   LOOP i4 = 1 to (numgroup-1). 
      COMPUTE m2int(:,i4) = mcent&*data(:,i4+1). 
   END LOOP. 
   COMPUTE datam2 = {MAKE(N,1,1), mcent, m2int}. 
   DO IF covtog = 1. 
      COMPUTE datam2 = {datam2, data(:,(2*numgroup+1):ncol(data))}.
   END IF. 
   COMPUTE yestm2 = datam2*inv(t(datam2)*datam2)*t(datam2)*y. 
   CORR var1 = y /var2 = yestm2.
   COMPUTE  ycorm2 = rPears. 
   COMPUTE redr2 = ycorm2**2.
!enddefine. 

define PROBE (min = !charend('/') /max = !charend('/')). 
COMPUTE jump = (!max - !min)/iter. 
COMPUTE dim = (transtog=0)*(iter+1) + (transtog=1)*(iter-1). 
COMPUTE tempres = MAKE(dim, 8, -999). 
DO IF transtog = 0. 
   COMPUTE i2 = !min. 
ELSE IF transtog = 1. 
   COMPUTE i2 = !min+jump. 
END IF.  
LOOP i = 1 to dim. 
   RedR center = i2. 
   COMPUTE tempres(i, 1:2) = {i2, redr2}.
   COMPUTE i2 = i2 + jump. 
END LOOP.  
COMPUTE tempres(:,3) = fullr2 - tempres(:,2). 
COMPUTE tempres(:,4) = (dffull*tempres(:,3))&/(dfred*(1-fullr2)).  
COMPUTE tempres(:,5) = 1-FCDF(tempres(:,4), dfred, dffull). 
!enddefine. 

define OGRS (vars = !charend('/') /x = !charend('/') /m = !charend('/') /y = !charend('/') 
           /conf = !charend('/') !default(95) /convcrit = !charend('/') !default(.00000001) /decimals = !charend('/') !default(F10.4) /iter = !charend('/') !default(0)). 
set mxloop = 100000000.
matrix. 
COMPUTE runnotes = MAKE(100,1,0). 
COMPUTE runcount = 1. 
COMPUTE criterr = 0. 
GET allvarsm /variables = !vars /names = allnames /missing = 999. 
GET allvars /variables = !vars /names = allnames /missing = OMIT. 
DO IF (nrow(allvars) <> nrow(allvarsm)). 
   COMPUTE runnotes(runcount, 1) = 1. 
   COMPUTE runcount = runcount +1. 
   COMPUTE missing = nrow(allvarsm) - nrow(allvars). 
END IF. 

GET xdat /variables = !x /names = xname /missing = OMIT. 
DO IF (ncol(xdat) > 1). 
   COMPUTE runnotes(runcount,1) = 2. 
   COMPUTE runcount = runcount +1. 
   COMPUTE criterr = 1. 
   COMPUTE xdat = xdat(:,1). 
   COMPUTE xname = xname(1,1). 
END IF. 
GET mdat /variables = !m /names = mname /missing = OMIT. 
DO IF (ncol(mdat) > 1). 
   COMPUTE runnotes(runcount,1) = 3. 
   COMPUTE runcount = runcount +1. 
   COMPUTE criterr = 1. 
   COMPUTE mdat = mdat(:,1). 
   COMPUTE mname = mname(1,1). 
END IF. 
GET ydat /variables = !y /names = yname /missing = OMIT. 
DO IF (ncol(ydat) > 1). 
   COMPUTE runnotes(runcount,1) = 4. 
   COMPUTE runcount = runcount +1. 
   COMPUTE criterr = 1. 
   COMPUTE ydat = ydat(:,1). 
   COMPUTE yname = yname(1,1). 
END IF. 

DO IF (criterr = 0). 
   LOOP i = 1 to ncol(allvars)-1. 
      LOOP j = i+1 to ncol(allvars). 
         DO IF (allnames(1,i) = allnames(1,j)). 
            COMPUTE runnotes(runcount,1) = 8. 
            COMPUTE runcount = runcount +1. 
            COMPUTE criterr = 1. 
         END IF. 
      END LOOP. 
   END LOOP. 
END IF. 

DO IF ((xname = mname) OR (xname = yname) OR (mname = yname)). 
   COMPUTE runnotes(runcount,1) = 8.  
   COMPUTE runcount = runcount +1. 
   COMPUTE criterr = 1. 
END IF. 

COMPUTE zero = MAKE(nrow(allvars),1,0).
LOOP i = 1 to (ncol(allvars)-1).
   LOOP j = i+1 to ncol(allvars). 
      COMPUTE diff = allvars(:,i) - allvars(:,j). 
      COMPUTE copy = (diff = zero).
      DO IF (all(copy) = 1). 
         COMPUTE copyname = {allnames(1,i), allnames(1,j)}. 
         COMPUTE runnotes(runcount,1) = 9. 
         COMPUTE runcount = runcount +1. 
         COMPUTE criterr = 1. 
      END IF. 
   END LOOP. 
END LOOP. 
   

COMPUTE convcrit = !convcrit. 
DO IF (convcrit < .00000001). 
   COMPUTE convcrit = .00000001. 
   COMPUTE runnotes(runcount, 1) = 6. 
   COMPUTE runcount = runcount +1. 
END IF. 

COMPUTE conf = !conf. 
DO IF ((conf < 50) OR (conf > 99.99999)). 
   COMPUTE conf = 95. 
   COMPUTE runnotes(runcount, 1) = 5. 
   COMPUTE runcount = runcount +1. 
END IF. 



DO IF (criterr <> 1). 

   COMPUTE alpha = 1-(conf/100). 
   COMPUTE covtog = (ncol(allvars) - 3 > 0). 
   COMPUTE N = nrow(allvars). 
   DO IF covtog =1. 
   COMPUTE covcount = 1. 
   COMPUTE cov = MAKE(N, ncol(allvars) - 3, 999).
   COMPUTE covname = MAKE(1, ncol(allvars)-3, 999). 
   END IF. 

   LOOP i = 1 to ncol(allnames). 
      DO IF (allnames(:,i) = xname). 
      COMPUTE xcol = i. 
      END IF. 
   END LOOP. 
   COMPUTE allvars(GRADE(allvars(:,xcol)),:) = allvars. 

   LOOP i = 1 to ncol(allnames). 
   DO IF (allnames(:,i) = xname). 
   COMPUTE x = allvars(:,i). 
   ELSE IF (allnames(:,i) = mname). 
   COMPUTE m = allvars(:,i). 
   ELSE IF (allnames(:,i) = yname). 
   COMPUTE y = allvars(:,i). 
   ELSE. 
   DO IF covtog = 1. 
   COMPUTE cov(:,covcount) = allvars(:,i). 
   COMPUTE covname(:,covcount) = allnames(:,i). 
   COMPUTE covcount = covcount +1. 
   END IF. 
   END IF. 
   END LOOP. 

   COMPUTE designX = design(x). 
   COMPUTE numgroup = ncol(designX). 
   COMPUTE designX = designX(:,1:(numgroup-1)). 
   COMPUTE xmat = MAKE(ncol(designX)+1, ncol(designX)+1, -999). 
   LOOP kloop = 1 to ncol(designX). 
   LOOP i = 1 to N. 
   DO IF (designx(i,kloop) = 1). 
   COMPUTE xmat(kloop,1) = x(i,1). 
   END IF. 
   END LOOP IF xmat(kloop,1) <> -999. 
   END LOOP. 
   LOOP i = 1 to N. 
   DO IF all(designx(i,:)=0). 
   COMPUTE xmat(ncol(designX)+1,:) = {x(i,1), MAKE(1,ncol(designX), 0)}. 
   END IF. 
   END LOOP if xmat(ncol(designX)+1,1) <> -999). 
   COMPUTE xmat(1:(numgroup-1),2:numgroup) = Ident(numgroup-1). 
   COMPUTE prodcol = MAKE(N, ncol(designX), 999). 
   LOOP i = 1 to ncol(designX). 
   COMPUTE prodcol(:,i) = designX(:,i)&*m. 
   END LOOP. 
   DO IF covtog = 0. 
   COMPUTE data = {MAKE(N,1,1), designX, m, prodcol}. 
   ELSE IF covtog = 1. 
   COMPUTE data = {MAKE(N,1,1), designX, m, prodcol, cov}. 
   END IF.
   DO IF (!iter <> 0). 
      COMPUTE iter = abs(rnd(!iter)). 
      DO IF (!iter <> iter). 
         COMPUTE runnotes(runcount,1) = 7. 
         COMPUTE runcount = runcount +1. 
      END IF. 
   ELSE IF (!iter = 0). 
      COMPUTE iter = 50+10*numgroup.  
   END IF. 
   COMPUTE yest = data*inv(t(data)*data)*t(data)*y. 
   CORR var1 = y /var2 = yest.
   COMPUTE  ycor = rPears. 
   COMPUTE fullr2 = ycor**2.
   COMPUTE dffull = N - ncol(data). 
   COMPUTE dfred = numgroup - 1. 
   COMPUTE Ffull = (fullr2*dffull)/((1-fullr2)*(ncol(data)-1)).
   COMPUTE pfull = 1-FCDF(Ffull, (ncol(data)-1), dffull). 
   COMPUTE modres = MAKE(ncol(data), 6, -999). 
   COMPUTE modres(:,1) = inv(t(data)*data)*t(data)*y. 
   COMPUTE ssr = csum((y - yest)&**2). 
   COMPUTE msr = ssr/(N-ncol(data)). 
   COMPUTE semat = (msr*inv(t(data)*data)).
   COMPUTE modres(:,2) = (diag(semat))&**(1/2). 
   COMPUTE modres(:,3) = modres(:,1)&/modres(:,2). 
   COMPUTE modres(:,4) = 2*(1-tcdf(abs(modres(:,3)),dffull)). 
   COMPUTE temp = alpha/2.
   CDFINVT p = temp /df = dffull. 
   COMPUTE tcrit = toutput. 
   COMPUTE modres(:,5) = modres(:,1) - tcrit*modres(:,2). 
   COMPUTE modres(:,6) = modres(:,1) + tcrit*modres(:,2). 
   DO IF covtog = 0. 
   COMPUTE dataint = {MAKE(N,1,1), designX, m}. 
   ELSE IF covtog = 1. 
   COMPUTE dataint = {MAKE(N,1,1), designX, m, cov}. 
   END IF. 
   COMPUTE yestint = dataint*inv(t(dataint)*dataint)*t(dataint)*y. 
   CORR var1 = y /var2 = yestint.
   COMPUTE  ycorint = rPears. 
   COMPUTE r2int = ycorint**2.
   COMPUTE rchint = fullr2 - r2int. 
   COMPUTE Fint = (dffull*rchint)&/(dfred*(1-fullr2)).  
   COMPUTE pint = 1-FCDF(Fint, dfred, dffull). 
   COMPUTE intres = {rchint, Fint, dfred, dffull, pint}. 
   COMPUTE transtog = 0. 
   COMPUTE minM = cmin(m). 
   COMPUTE maxM = cmax(m). 
   PROBE min = minM /max = maxM. 
   COMPUTE results = tempres.
   COMPUTE OGres = tempres.  
   COMPUTE results(nrow(results),6:7) = {0,0}. 
   COMPUTE i3 = 1. 
   LOOP IF i3 <= nrow(results). 
   DO IF i3 < nrow(results). 
   COMPUTE results(i3, 6) = 1*(results(i3,4) < results(i3+1,4)) - 1*(results(i3,4) > results(i3+1,4)).
   COMPUTE results(i3, 7) = -1*((results(i3,5) < alpha) AND (results(i3+1,5) > alpha)) + 1*((results(i3,5) > alpha) AND (results(i3+1, 5) < alpha)). 
   END IF. 
   COMPUTE results(i3,8) = (abs(results(i3,5) - alpha) < convcrit). 
   DO IF i3 = nrow(results). 
   COMPUTE transcnv = 0. 
   ELSE IF i3 = 1. 
   COMPUTE transcnv = ((results(i3,7) = 1)AND((results(i3,8)=1)OR(abs(results(i3+1,5) - alpha) < convcrit))). 
   ELSE.    
   COMPUTE trnscnv1 = ((results(i3,7) = 1) AND ((results(i3,8)=1) OR (abs(results(i3+1,5) - alpha) < convcrit))).
   COMPUTE trnscnv2 = ((results(i3,7) = -1) AND ((results(i3,8)=1) OR (abs(results(i3+1,5) - alpha) < convcrit))). 
   COMPUTE transcnv = (trnscnv1 = 1) OR (trnscnv2 = 1). 
   END IF. 
   DO IF ((abs(results(i3,7))=1)AND(transcnv = 0)). 
   COMPUTE trnsindx = i3. 
   COMPUTE transtog = 1. 
   COMPUTE minmtran = mmin({results(i3+1,1), results(i3,1)}). 
   COMPUTE maxmtran = mmax({results(i3+1,1), results(i3,1)}). 
   PROBE min = minmtran /max = maxmtran. 
   COMPUTE results = {results; tempres}. 
   COMPUTE results(GRADE(results(:,1)),:) = results. 
   ELSE. 
   COMPUTE i3 = i3+1. 
   END IF. 
   END LOOP. 
   COMPUTE numJN = 1*(results(nrow(results),8) =1) + 1*((results(1,8) = 1) AND (results(1,7) <> 1)) + csum(abs(results(:,7))). 
   DO IF numJN > 0. 
   COMPUTE JNSoln = MAKE(numJN,1, -999). 
   COMPTUE JNIndx = MAKE(numJN, 1, -999). 
   COMPUTE slncnt = 1. 
   DO IF results(nrow(results),8) = 1. 
   COMPUTE JNSoln(1,1) = results(nrow(results),1). 
   COMPUTE JNIndx(1,1) = nrow(results). 
   COMPUTE slncnt = slncnt +1. 
   END IF. 
   LOOP i1 = 1 to nrow(results). 
   DO IF abs(results(i1,7)) = 1.
   COMPUTE abvblw = {results(i1,1), abs(results(i1,5)-alpha); results(i1+1, 1), abs(results(i1+1,5) - alpha)}. 
   COMPUTE minval = GRADE(abvblw(:,2)). 
   COMPUTE indxtog = all(abvblw(GRADE(abvblw(:,2)),:) = abvblw). 
   DO IF (indxtog = 1). 
   COMPUTE JNIndx(slncnt,1) = i1. 
   ELSE. 
   COMPUTE JNIndx(slncnt,1) = i1+1. 
   END IF. 
   COMPUTE abvblw(GRADE(abvblw(:,2)),:) = abvblw.    
   COMPUTE JNSoln(slncnt,1) = abvblw(1,1).
   COMPUTE slncnt = slncnt+1. 
   END IF. 
   END LOOP. 
   END IF. 
END IF. 
PRINT /title = "******************* OGRS Procedure for SPSS Version 1.2 ********************".
PRINT /title = "                           Written by Amanda Montoya       ".
PRINT /title = "                       Documentation available by request ".
PRINT /title = "*****************************************************************************".
DO IF (criterr <>1). 
   COMPUTE varrlabs = {'X =', 'M = ', 'Y = '}. 
   PRINT {xname; mname; yname} /title = "Variables:" /rnames = varrlabs /format = A8. 
   DO IF covtog = 1. 
   PRINT {covname} /title = "Statistical Controls:" /format = A8. 
   END IF. 
   COMPUTE dummylab = {"D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"}. 
   COMPUTE xmatlab = {xname, dummylab(1:(numgroup-1))}.
   PRINT xmat /title = "Dummy Variable Coding Scheme:" /cnames = xmatlab. 
   PRINT N /title = "Sample size:".
   PRINT {yname} /title = "******************************************************************************" /rlabels = "Outcome:" /format = A8.
   COMPUTE modsum = {sqrt(fullr2), fullr2, Ffull, (ncol(data)-1), dffull, pfull}. 
   PRINT modsum /title = "Model Summary" /clabels = "R", "R-sq", "F" , "df1" , "df2", "p" /format = !decimals.
   COMPUTE intlab = {"Int1", "Int2", "Int3", "Int4", "Int5", "Int6", "Int7", "Int8", "Int9"}. 
   COMPUTE modlabs = {"constant", dummylab(1,1:(numgroup-1)), mname, intlab(1,1:(numgroup-1))}. 
   DO IF (covtog = 1). 
   COMPUTE modlabs = {modlabs, covname}. 
   END IF.  
   PRINT modres /title "Model" /rnames = modlabs /clabels = "coeff" , "SE", "t", "p", "LLCI", "ULCI" /format = !decimals.
   COMPUTE intmat = MAKE((numgroup-1), 5, -999). 
   COMPUTE intmat(:,1) = t(intlab(1, 1:(numgroup-1))).
   COMPUTE intmat(:,2) = MAKE((numgroup-1), 1, "=").  
   COMPUTE intmat(:,3) = t(dummylab(1, 1:(numgroup-1))).
   COMPUTE intmat(:,4) = MAKE((numgroup-1), 1, "X"). 
   COMPUTE intmat(:,5) = MAKE((numgroup-1), 1, mname). 
   PRINT intmat /title = "Interactions:" /format = A8. 
   PRINT intres /title = "R-square increase due to interaction(s):" /clabels = "R2-chng" "F" "df1" "df2" "p" /format = !decimals. 
   PRINT /title = "************************* JOHNSON-NEYMAN TECHNIQUE **************************".
   DO IF (iter > 10). 
   COMPUTE last = nrow(OGres).
   COMPUTE rjump = rnd(last/20).
   COMPTUE rowsel = 1.
   COMPUTE rcount = 1+rjump.
   LOOP IF (rcount <= last).
   COMPUTE rowsel = {rowsel, rcount}. 
   COMPUTE rcount = rcount + rjump. 
   END LOOP.
   DO IF (rcount - rjump <> last). 
   COMPUTE rowsel = {rowsel, last}. 
   END IF. 
   END IF. 
   COMPUTE JNtabnam = {mname, "R2-chng", "F", "p"}. 
   DO IF numJN > 0. 
   PRINT JNSoln /title = "Moderator value(s) defining Johnson-Neyman boundaries of significance:" /format = !decimals.
   DO IF (iter > 10). 
   COMPUTE JNouttab = {OGres(rowsel,:); results(JNIndx, :)}.
   ELSE. 
   COMPUTE JNouttab = {OGres(:,:);results(JNIndx,:)}. 
   END IF. 
   COMPUTE JNouttab(GRADE(JNouttab(:,1)),:) = JNouttab.
   COMPUTE JNouttab = JNouttab(:,{1,3,4,5}). 
   PRINT JNouttab /title = "Conditional effect of X on Y at values of the moderator (M)" /cnames = JNtabnam /format = !decimals. 
   ELSE. 
   PRINT /title = "No Johnson-Neyman bounds found within range of observed data".
   DO IF (iter > 10). 
   COMPUTE JNouttab = OGres(rowsel,{1,3,4,5}). 
   ELSE. 
   COMPUTE JNouttab = OGres(:,{1,3,4,5}). 
   END IF. 
   PRINT JNouttab /title = "Conditional effect of X on Y at values of the moderator (M)" /cnames = JNtabnam /format = !decimals. 
   END IF. 
END IF. 
PRINT /title = "************************ ANALYSIS NOTES AND WARNINGS ************************".
DO IF (criterr <>1). 
   PRINT nrow(results) /title = "Number of points probed in Johnson-Neyman algorithm".
END IF. 
LOOP i = 1 to nrow(runnotes). 
   DO IF (runnotes(i,1) = 1). 
      print missing /title = "NOTE: Some cases were deleted due to missing data. The number of cases was:".
   ELSE IF (runnotes(i,1) = 2). 
      print /title = "ERROR: Only one X variable can be specified in the X = list.".
   ELSE IF (runnotes(i,1) = 3). 
      print /title = "ERROR: Only one M variable can be specified in the M = list.". 
   ELSE IF (runnotes(i,1) = 4). 
      print /title = "ERROR: Only one Y variable can be specified in the Y = list.".
   ELSE IF (runnotes(i,1) = 5). 
      print /title = "NOTE: The confidence level specified was not between 50 and 99.9999.".
      print /title = "      Level of confidence was adjusted to 95%.".
   ELSE IF (runnotes(i,1) = 6). 
      print /title = "NOTE: Convergence criteria specified was not greater than or equal to .00000001.". 
      print /title = "      Level of convergence criteria was adjusted to .00000001.".
   ELSE IF (runnotes(i,1) = 7). 
      print /title = "NOTE: An invalid number of initial iterations was specified.".
      print iter /title = "      The number of initial iterations used was:". 
   ELSE IF (runnotes(i,1) = 8). 
      print /title = "ERROR: All variables must be unique. No variables can be the same among X, M, and Y". 
   ELSE IF (runnotes(i,1) = 9). 
      print copyname /title = "ERROR: Two of the specified variables are copies. The variables names are:" /format = A8. 
   END IF. 
END LOOP. 
PRINT /title = "*****************************************************************************".
end matrix. 
!enddefine. 
restore. 


