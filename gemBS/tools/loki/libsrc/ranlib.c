#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ranlib.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

/************************************************************************
FTNSTOP:
Prints msg to standard error and then exits
************************************************************************/
static void ftnstop(char* msg)
{
  if (msg != NULL) (void)fprintf(stderr,"%s\n",msg);
  exit(EXIT_SUCCESS);
}




static double fsign( const double num, const double sign)
{
	return ((sign>DBL_EPSILON && num<DBL_EPSILON)||(sign<DBL_EPSILON && num>DBL_EPSILON))?-num:num;
}

double genbet(const double aa,const double bb)
/**********************************************************************
     double genbet(double aa,double bb)
               GeNerate BETa random deviate
                              Function
     Returns a single random deviate from the beta distribution with
     parameters A and B.  The density of the beta is
               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
                              Arguments
     aa --> First parameter of the beta distribution
       
     bb --> Second parameter of the beta distribution
       
                              Method
     R. C. H. Cheng
     Generating Beta Variatew with Nonintegral Shape Parameters
     Communications of the ACM, 21:317-322  (1978)
     (Algorithms BB and BC)
**********************************************************************
*/
{
/* JJV changed expmax (log(1.0E38)==87.49823), and added minlog */
#define expmax 87.49823
#define infnty 1.0E38
#define minlog 1.0E-37
static double olda = -1.0E37;
static double oldb = -1.0E37;
static double gbet,a,alpha,b,beta,delta,gamm,k1,k2,r,s,t,u1,u2,v,w,y,z;
static int qsame;

    qsame = fabs(olda-aa)<DBL_EPSILON && fabs(oldb-bb)<DBL_EPSILON;
    if(qsame) goto S20;
    if(!(aa < minlog || bb < minlog)) goto S10;
    (void)fputs(" AA or BB < 1.0E-37 in GENBET - Abort!\n",stderr);
    (void)fprintf(stderr," AA: %16.6E BB %16.6E\n",aa,bb);
    exit(EXIT_FAILURE);
S10:
    olda = aa;
    oldb = bb;
S20:
    if(!(min(aa,bb) > 1.0)) goto S100;
/*
     Alborithm BB
     Initialize
*/
    if(qsame) goto S30;
    a = min(aa,bb);
    b = max(aa,bb);
    alpha = a+b;
    beta = sqrt((alpha-2.0)/(2.0*a*b-alpha));
    gamm = a+1.0/beta;
S30:
S40:
    u1 = safe_ranf();
/*
     Step 1
*/
    u2 = safe_ranf();
    v = beta*log(u1/(1.0-u1));
/* JJV altered this */
    if(v > expmax) goto S55;
/*
 * JJV added checker to see if a*exp(v) will overflow
 * JJV S50 _was_ w = a*exp(v); also note here a > 1.0
 */
    w = exp(v);
    if(w > infnty/a) goto S55;
    w *= a;
    goto S60;
S55:
    w = infnty;
S60:
    z = pow(u1,2.0)*u2;
    r = gamm*v-1.3862944;
    s = a+r-w;
/*
     Step 2
*/
    if(s+2.609438 >= 5.0*z) goto S70;
/*
     Step 3
*/
    t = log(z);
    if(s > t) goto S70;
/*
 *   Step 4
 *
 *    JJV added checker to see if log(alpha/(b+w)) will 
 *    JJV overflow.  If so, we count the log as -INF, and
 *    JJV consequently evaluate conditional as true, i.e.
 *    JJV the algorithm rejects the trial and starts over
 *    JJV May not need this here since alpha > 2.0
 */
    if(alpha/(b+w) < minlog) goto S40;
    if(r+alpha*log(alpha/(b+w)) < t) goto S40;
S70:
/*
     Step 5
*/
    if(!(aa == a)) goto S80;
    gbet = w/(b+w);
    goto S90;
S80:
    gbet = b/(b+w);
S90:
    goto S230;
S100:
/*
     Algorithm BC
     Initialize
*/
    if(qsame) goto S110;
    a = max(aa,bb);
    b = min(aa,bb);
    alpha = a+b;
    beta = 1.0/b;
    delta = 1.0+a-b;
    k1 = delta*(1.38889E-2+4.16667E-2*b)/(a*beta-0.777778);
    k2 = 0.25+(0.5+0.25/delta)*b;
S110:
S120:
    u1 = safe_ranf();
/*
     Step 1
*/
    u2 = safe_ranf();
    if(u1 >= 0.5) goto S130;
/*
     Step 2
*/
    y = u1*u2;
    z = u1*y;
    if(0.25*u2+z-y >= k1) goto S120;
    goto S170;
S130:
/*
     Step 3
*/
    z = pow(u1,2.0)*u2;
    if(!(z <= 0.25)) goto S160;
    v = beta*log(u1/(1.0-u1));
/*
 *    JJV instead of checking v > expmax at top, I will check
 *    JJV if a < 1, then check the appropriate values
 */
    if(a > 1.0) goto S135;
/*   JJV a < 1 so it can help out if exp(v) would overflow */
    if(v > expmax) goto S132;
    w = a*exp(v);
    goto S200;
S132:
    w = v + log(a);
    if(w > expmax) goto S140;
    w = exp(w);
    goto S200;
S135:
/*   JJV in this case a > 1 */
    if(v > expmax) goto S140;
    w = exp(v);
    if(w > infnty/a) goto S140;
    w *= a;
    goto S200;
S140:
    w = infnty;
    goto S200;
/*
 * JJV old code
 *    if(!(v > expmax)) goto S140;
 *    w = infnty;
 *    goto S150;
 *S140:
 *    w = a*exp(v);
 *S150:
 *    goto S200;
 */
S160:
    if(z >= k2) goto S120;
S170:
/*
     Step 4
     Step 5
*/
    v = beta*log(u1/(1.0-u1));
/*   JJV same kind of checking as above */
    if(a > 1.0) goto S175;
/* JJV a < 1 so it can help out if exp(v) would overflow */
    if(v > expmax) goto S172;
    w = a*exp(v);
    goto S190;
S172:
    w = v + log(a);
    if(w > expmax) goto S180;
    w = exp(w);
    goto S190;
S175:
/* JJV in this case a > 1.0 */
    if(v > expmax) goto S180;
    w = exp(v);
    if(w > infnty/a) goto S180;
    w *= a;
    goto S190;
S180:
    w = infnty;
/*
 *   JJV old code
 *    if(!(v > expmax)) goto S180;
 *    w = infnty;
 *    goto S190;
 *S180:
 *    w = a*exp(v);
 */
S190:
/*
 * JJV here we also check to see if log overlows; if so, we treat it
 * JJV as -INF, which means condition is true, i.e. restart
 */
    if(alpha/(b+w) < minlog) goto S120;
    if(alpha*(log(alpha/(b+w))+v)-1.3862944 < log(z)) goto S120;
S200:
/*
     Step 6
*/
    if(!(a == aa)) goto S210;
    gbet = w/(b+w);
    goto S220;
S210:
    gbet = b/(b+w);
S230:
S220:
    return gbet;
#undef expmax
#undef infnty
#undef minlog
}
double genchi(const double df)
/**********************************************************************
     double genchi(double df)
                Generate random value of CHIsquare variable
                              Function
     Generates random deviate from the distribution of a chisquare
     with DF degrees of freedom random variable.
                              Arguments
     df --> Degrees of freedom of the chisquare
            (Must be positive)
       
                              Method
     Uses relation between chisquare and gamma.
**********************************************************************
*/
{
static double gchi;

    if(!(df <= 0.0)) goto S10;
    (void)fputs(" DF <= 0 in GENCHI - ABORT\n",stderr);
    (void)fprintf(stderr," Value of DF: %16.6E\n",df);
    exit(EXIT_FAILURE);
S10:
/*
 * JJV changed the code to call SGAMMA directly
 *    genchi = 2.0*gengam(1.0,df/2.0); <- OLD
 */
    gchi = 2.0*sgamma(df/2.0);
    return gchi;
}

double genexp(double const av)
/*
**********************************************************************
     double genexp(double av)
                    GENerate EXPonential random deviate
                              Function
     Generates a single random deviate from an exponential
     distribution with mean AV.
                              Arguments
     av --> The mean of the exponential distribution from which
            a random deviate is to be generated.
        JJV (av >= 0)
                              Method
     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling From the
               Exponential and Normal Distributions.
               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
**********************************************************************
*/
{
static double gexp;

/* JJV added check that av >= 0 */
    if(av >= 0.0) goto S10;
    (void)fputs(" AV < 0 in GENEXP - ABORT\n",stderr);
    (void)fprintf(stderr," Value of AV: %16.6E\n",av);
    exit(EXIT_FAILURE);
S10:
    gexp = sexpo()*av;
    return gexp;
}
double genf(const double dfn,const double dfd)
/*
**********************************************************************
     double genf(double dfn,double dfd)
                GENerate random deviate from the F distribution
                              Function
     Generates a random deviate from the F (variance ratio)
     distribution with DFN degrees of freedom in the numerator
     and DFD degrees of freedom in the denominator.
                              Arguments
     dfn --> Numerator degrees of freedom
             (Must be positive)
     dfd --> Denominator degrees of freedom
             (Must be positive)
                              Method
     Directly generates ratio of chisquare variates
**********************************************************************
*/
{
static double gf,xden,xnum;

    if(!(dfn <= 0.0 || dfd <= 0.0)) goto S10;
    (void)fputs(" Degrees of freedom nonpositive in GENF - abort!\n",stderr);
    (void)fprintf(stderr," DFN value: %16.6E DFD value: %16.6E\n",dfn,dfd);
    exit(EXIT_FAILURE);
S10:
/*
 * JJV changed this to call SGAMMA directly
 *
 *     GENF = ( GENCHI( DFN ) / DFN ) / ( GENCHI( DFD ) / DFD )
 *   xnum = genchi(dfn)/dfn; <- OLD
 *   xden = genchi(dfd)/dfd; <- OLD
 */
    xnum = 2.0*sgamma(dfn/2.0)/dfn;
    xden = 2.0*sgamma(dfd/2.0)/dfd;
/*
 * JJV changed constant to prevent underflow at compile time.
 *   if(!(xden <= 9.999999999998E-39*xnum)) goto S20;
 */
    if(!(xden <= 1.0E-37*xnum)) goto S20;
    (void)fputs(" GENF - generated numbers would cause overflow\n",stderr);
    (void)fprintf(stderr," Numerator %16.6E Denominator %16.6E\n",xnum,xden);
/*
 * JJV changed next 2 lines to reflect constant change above in the
 * JJV truncated value returned.
 *   (void)fputs(" GENF returning 1.0E38\n",stderr);
 *   genf = 1.0E38;
 */
    (void)fputs(" GENF returning 1.0E37\n",stderr);
    gf = 1.0E37;
    goto S30;
S20:
    gf = xnum/xden;
S30:
    return gf;
}
double gengam(const double a,const double r)
/*
**********************************************************************
     double gengam(double a,double r)
           GENerates random deviates from GAMma distribution
                              Function
     Generates random deviates from the gamma distribution whose
     density is
          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
                              Arguments
     a --> Location parameter of Gamma distribution
     JJV   (a > 0)
     r --> Shape parameter of Gamma distribution
     JJV   (r > 0)
                              Method
     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               (Case R >= 1.0)
               Ahrens, J.H. and Dieter, U.
               Generating Gamma Variates by a
               Modified Rejection Technique.
               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
     Algorithm GD
     JJV altered following to reflect argument ranges
               (Case 0.0 < R < 1.0)
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling from Gamma,
               Beta, Poisson and Binomial Distributions.
               Computing, 12 (1974), 223-246/
     Adapted algorithm GS.
**********************************************************************
*/
{
static double ggam;
/* JJV added argument checker */
    if(a > 0.0 && r > 0.0) goto S10;
    (void)fputs(" A or R nonpositive in GENGAM - abort!\n",stderr);
    (void)fprintf(stderr," A value: %16.6E R value: %16.6E\n",a,r);
    exit(EXIT_FAILURE);
S10:
    ggam = sgamma(r);
    ggam /= a;
    return ggam;
}

void genmul(const int n,const double *p,const int ncat,int *ix)
/*
**********************************************************************
 
     void genmul(int n,double *p,int ncat,int *ix)
     GENerate an observation from the MULtinomial distribution
                              Arguments
     N --> Number of events that will be classified into one of
           the categories 1..NCAT
     P --> Vector of probabilities.  P(i) is the probability that
           an event will be classified into category i.  Thus, P(i)
           must be [0,1]. Only the first NCAT-1 P(i) must be defined
           since P(NCAT) is 1.0 minus the sum of the first
           NCAT-1 P(i).
     NCAT --> Number of categories.  Length of P and IX.
     IX <-- Observation from multinomial distribution.  All IX(i)
            will be nonnegative and their sum will be N.
                              Method
     Algorithm from page 559 of
 
     Devroye, Luc
 
     Non-Uniform Random Variate Generation.  Springer-Verlag,
     New York, 1986.
 
**********************************************************************
*/
{
static double prob,ptot,sum;
static int i,icat,ntot;
    if(n < 0) ftnstop("N < 0 in GENMUL");
    if(ncat <= 1) ftnstop("NCAT <= 1 in GENMUL");
    ptot = 0.0F;
    for(i=0; i<ncat-1; i++) {
        if(*(p+i) < 0.0F) ftnstop("Some P(i) < 0 in GENMUL");
        if(*(p+i) > 1.0F) ftnstop("Some P(i) > 1 in GENMUL");
        ptot += *(p+i);
    }
    if(ptot > 0.99999F) ftnstop("Sum of P(i) > 1 in GENMUL");
/*
     Initialize variables
*/
    ntot = n;
    sum = 1.0F;
    for(i=0; i<ncat; i++) ix[i] = 0;
/*
     Generate the observation
*/
    for(icat=0; icat<ncat-1; icat++) {
        prob = *(p+icat)/sum;
        *(ix+icat) = ignbin(ntot,prob);
        ntot -= *(ix+icat);
	if(ntot <= 0) return;
        sum -= *(p+icat);
    }
    *(ix+ncat-1) = ntot;
/*
     Finished
*/
    return;
}

double gennch(const double df,const double xnonc)
/*
**********************************************************************
     double gennch(double df,double xnonc)
           Generate random value of Noncentral CHIsquare variable
                              Function
     Generates random deviate  from the  distribution  of a  noncentral
     chisquare with DF degrees  of freedom and noncentrality  parameter
     xnonc.
                              Arguments
     df --> Degrees of freedom of the chisquare
            (Must be >= 1.0)
     xnonc --> Noncentrality parameter of the chisquare
               (Must be >= 0.0)
                              Method
     Uses fact that  noncentral chisquare  is  the  sum of a  chisquare
     deviate with DF-1  degrees of freedom plus the  square of a normal
     deviate with mean XNONC and standard deviation 1.
**********************************************************************
*/
{
static double g;

    if(!(df < 1.0 || xnonc < 0.0)) goto S10;
    (void)fputs("DF < 1 or XNONC < 0 in GENNCH - ABORT\n",stderr);
    (void)fprintf(stderr,"Value of DF: %16.6E Value of XNONC: %16.6E\n",df,xnonc);
    exit(EXIT_FAILURE);
/* JJV changed code to call SGAMMA, SNORM directly */
S10:
    if(df >= 1.000001) goto S20;
/*
 * JJV case df == 1.0
 * gennch = pow(gennor(sqrt(xnonc),1.0),2.0); <- OLD
 */
    g = pow(snorm()+sqrt(xnonc),2.0);
    goto S30;
S20:
/*
 * JJV case df > 1.0
 * gennch = genchi(df-1.0)+pow(gennor(sqrt(xnonc),1.0),2.0); <- OLD
 */
    g = 2.0*sgamma((df-1.0)/2.0)+pow(snorm()+sqrt(xnonc),2.0);
S30:
    return g;
}

double gennf(const double dfn,const double dfd,const double xnonc)
/*
**********************************************************************
     double gennf(double dfn,double dfd,double xnonc)
           GENerate random deviate from the Noncentral F distribution
                              Function
     Generates a random deviate from the  noncentral F (variance ratio)
     distribution with DFN degrees of freedom in the numerator, and DFD
     degrees of freedom in the denominator, and noncentrality parameter
     XNONC.
                              Arguments
     dfn --> Numerator degrees of freedom
             (Must be >= 1.0)
     dfd --> Denominator degrees of freedom
             (Must be positive)
     xnonc --> Noncentrality parameter
               (Must be nonnegative)
                              Method
     Directly generates ratio of noncentral numerator chisquare variate
     to central denominator chisquare variate.
**********************************************************************
*/
{
static double gnf,xden,xnum;
static int qcond;

    /* JJV changed qcond, error message to allow dfn == 1.0 */
    qcond = dfn < 1.0 || dfd <= 0.0 || xnonc < 0.0;
    if(!qcond) goto S10;
    (void)fputs("In GENNF - Either (1) Numerator DF < 1.0 or\n",stderr);
    (void)fputs(" (2) Denominator DF <= 0.0 or\n",stderr);
    (void)fputs(" (3) Noncentrality parameter < 0.0\n",stderr);
    (void)fprintf(stderr,
      "DFN value: %16.6E DFD value: %16.6E XNONC value: \n%16.6E\n",dfn,dfd,
      xnonc);
    exit(EXIT_FAILURE);
S10:
/*
 * JJV changed the code to call SGAMMA and SNORM directly
 * GENNF = ( GENNCH( DFN, XNONC ) / DFN ) / ( GENCHI( DFD ) / DFD )
 * xnum = gennch(dfn,xnonc)/dfn; <- OLD
 * xden = genchi(dfd)/dfd; <- OLD
 */
    if(dfn >= 1.000001) goto S20;
/* JJV case dfn == 1.0, dfn is counted as exactly 1.0 */
    xnum = pow(snorm()+sqrt(xnonc),2.0);
    goto S30;
S20:
/* JJV case df > 1.0 */
    xnum = (2.0*sgamma((dfn-1.0)/2.0)+pow(snorm()+sqrt(xnonc),2.0))/dfn;
S30:
    xden = 2.0*sgamma(dfd/2.0)/dfd;
/*
 * JJV changed constant to prevent underflow at compile time.
 *   if(!(xden <= 9.999999999998E-39*xnum)) goto S40;
 */
    if(!(xden <= 1.0E-37*xnum)) goto S40;
    (void)fputs(" GENNF - generated numbers would cause overflow\n",stderr);
    (void)fprintf(stderr," Numerator %16.6E Denominator %16.6E\n",xnum,xden);
/*
 * JJV changed next 2 lines to reflect constant change above in the
 * JJV truncated value returned.
 *   (void)fputs(" GENNF returning 1.0E38\n",stderr);
 *   gnf = 1.0E38;
 */
    (void)fputs(" GENNF returning 1.0E37\n",stderr);
    gnf = 1.0E37;
    goto S50;
S40:
    gnf = xnum/xden;
S50:
    return gnf;
}
double gennor(const double av,const double sd)
/*
**********************************************************************
     double gennor(double av,double sd)
         GENerate random deviate from a NORmal distribution
                              Function
     Generates a single random deviate from a normal distribution
     with mean, AV, and standard deviation, SD.
                              Arguments
     av --> Mean of the normal distribution.
     sd --> Standard deviation of the normal distribution.
     JJV    (sd >= 0)
                              Method
     Renames SNORM from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               Ahrens, J.H. and Dieter, U.
               Extensions of Forsythe's Method for Random
               Sampling from the Normal Distribution.
               Math. Comput., 27,124 (Oct. 1973), 927 - 937.
**********************************************************************
*/
{
/* JJV added argument checker */
    if(sd >= 0.0) goto S10;
    (void)fputs(" SD < 0 in GENNOR - ABORT\n",stderr);
    (void)fprintf(stderr," Value of SD: %16.6E\n",sd);
    exit(EXIT_FAILURE);
S10:
    return sd*snorm()+av;
}
void genprm(int *iarray,const int larray)
/*
**********************************************************************
    void genprm(int *iarray,int larray)
               GENerate random PeRMutation of iarray
                              Arguments
     iarray <--> On output IARRAY is a random permutation of its
                 value on input
     larray <--> Length of IARRAY
**********************************************************************
*/
{
static int i,itmp,iwhich,D1,D2;

    for(i=1,D1=1,D2=(larray-i+D1)/D1; D2>0; D2--,i+=D1) {
        iwhich = ignuin(i,larray);
        itmp = *(iarray+iwhich-1);
        *(iarray+iwhich-1) = *(iarray+i-1);
        *(iarray+i-1) = itmp;
    }
}
double genunf(const double low,const double high)
/*
**********************************************************************
     double genunf(double low,double high)
               GeNerate Uniform Real between LOW and HIGH
                              Function
     Generates a real uniformly distributed between LOW and HIGH.
                              Arguments
     low --> Low bound (exclusive) on real value to be generated
     high --> High bound (exclusive) on real value to be generated
**********************************************************************
*/
{
    if(!(low > high)) goto S10;
    (void)fprintf(stderr,"LOW > HIGH in GENUNF: LOW %16.6E HIGH: %16.6E\n",low,high);
    (void)fputs("Abort\n",stderr);
    exit(EXIT_FAILURE);
S10:
    return low+(high-low)*safe_ranf();
}

int ignbin(const int n,const double pp)
/*
**********************************************************************
     int ignbin(int n,double pp)
                    GENerate BINomial random deviate
                              Function
     Generates a single random deviate from a binomial
     distribution whose number of trials is N and whose
     probability of an event in each trial is P.
                              Arguments
     n  --> The number of trials in the binomial distribution
            from which a random deviate is to be generated.
	    JJV (N >= 0)
     pp --> The probability of an event in each trial of the
            binomial distribution from which a random deviate
            is to be generated.
	    JJV (0.0 <= PP <= 1.0)
     ignbin <-- A random deviate yielding the number of events
                from N independent trials, each of which has
                a probability of event P.
                              Method
     This is algorithm BTPE from:
         Kachitvichyanukul, V. and Schmeiser, B. W.
         Binomial Random Variate Generation.
         Communications of the ACM, 31, 2
         (February, 1988) 216.
**********************************************************************
     SUBROUTINE BTPEC(N,PP,ISEED,JX)
     BINOMIAL RANDOM VARIATE GENERATOR
     MEAN .LT. 30 -- INVERSE CDF
       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
       USABLE ALGORITHM.
     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
       "BINOMIAL RANDOM VARIATE GENERATION,"
       COMMUNICATIONS OF THE ACM, FORTHCOMING
     WRITTEN:  SEPTEMBER 1980.
       LAST REVISED:  MAY 1985, JULY 1987
     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
                           GENERATOR
     ARGUMENTS
       N : NUMBER OF BERNOULLI TRIALS            (INPUT)
       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
     VARIABLES
       PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
       NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
       M:  INTEGER VALUE OF THE CURRENT MODE
       FM:  doubleING POINT VALUE OF THE CURRENT MODE
       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
       P1:  AREA OF THE TRIANGLE
       C:  HEIGHT OF THE PARALLELOGRAMS
       XM:  CENTER OF THE TRIANGLE
       XL:  LEFT END OF THE TRIANGLE
       XR:  RIGHT END OF THE TRIANGLE
       AL:  TEMPORARY VARIABLE
       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
       P2:  AREA OF THE PARALLELOGRAMS
       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
           FROM THE REGION
       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
           REJECT THE CANDIDATE VALUE
       IX:  INTEGER CANDIDATE VALUE
       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
           AND A doubleING POINT IX IN THE ACCEPT/REJECT LOGIC
       K:  ABSOLUTE VALUE OF (IX-M)
       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
           ALSO USED IN THE INVERSE TRANSFORMATION
       R: THE RATIO P/Q
       G: CONSTANT USED IN CALCULATION OF PROBABILITY
       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
            OF F WHEN IX IS GREATER THAN M
       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
             CALCULATION OF F WHEN IX IS LESS THAN M
       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
       YNORM: LOGARITHM OF NORMAL BOUND
       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
       USED IN THE FINAL ACCEPT/REJECT TEST
       QN: PROBABILITY OF NO SUCCESS IN N TRIALS
     REMARK
       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
       ARE NOT INVOLVED.
     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/
{
/* JJV changed initial values to ridiculous values */
static double psave = -1.0E37;
static int nsave = -214748365;
static int i,ix,ix1,k,m,mp,T1;
static double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
    x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

    if(pp != psave) goto S10;
    if(n != nsave) goto S20;
    if(xnp < 30.0) goto S150;
    goto S30;
S10:
/*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
JJV added checks to ensure 0.0 <= PP <= 1.0
*/
    if(pp < 0.0F) ftnstop("PP < 0.0 in IGNBIN");
    if(pp > 1.0F) ftnstop("PP > 1.0 in IGNBIN");
    psave = pp;
    p = min(psave,1.0-psave);
    q = 1.0-p;
S20:
/*
JJV added check to ensure N >= 0
*/
    if(n < 0L) ftnstop("N < 0 in IGNBIN");
    xnp = n*p;
    nsave = n;
    if(xnp < 30.0) goto S140;
    ffm = xnp+p;
    m = (int)ffm;
    fm = (double)m;
    xnpq = xnp*q;
    p1 = (int) (2.195*sqrt(xnpq)-4.6*q)+0.5;
    xm = fm+0.5;
    xl = xm-p1;
    xr = xm+p1;
    c = 0.134+20.5/(15.3+fm);
    al = (ffm-xl)/(ffm-xl*p);
    xll = al*(1.0+0.5*al);
    al = (xr-ffm)/(xr*q);
    xlr = al*(1.0+0.5*al);
    p2 = p1*(1.0+c+c);
    p3 = p2+c/xll;
    p4 = p3+c/xlr;
S30:
/*
*****GENERATE VARIATE
*/
    u = safe_ranf()*p4;
    v = safe_ranf();
/*
     TRIANGULAR REGION
*/
    if(u > p1) goto S40;
    ix = (int)(xm-p1*v+u);
    goto S170;
S40:
/*
     PARALLELOGRAM REGION
*/
    if(u > p2) goto S50;
    x = xl+(u-p1)/c;
    v = v*c+1.0-ABS(xm-x)/p1;
    if(v > 1.0 || v <= 0.0) goto S30;
    ix = (int)x;
    goto S70;
S50:
/*
     LEFT TAIL
*/
    if(u > p3) goto S60;
    ix = (int)(xl+log(v)/xll);
    if(ix < 0) goto S30;
    v *= ((u-p2)*xll);
    goto S70;
S60:
/*
     RIGHT TAIL
*/
    ix = (int)(xr-log(v)/xlr);
    if(ix > n) goto S30;
    v *= ((u-p3)*xlr);
S70:
/*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
    k = ABS(ix-m);
    if(k > 20 && k < xnpq/2-1) goto S130;
/*
     EXPLICIT EVALUATION
*/
    f = 1.0;
    r = p/q;
    g = (n+1)*r;
    T1 = m-ix;
    if(T1 < 0) goto S80;
    else if(T1 == 0) goto S120;
    else  goto S100;
S80:
    mp = m+1;
    for(i=mp; i<=ix; i++) f *= (g/i-r);
    goto S120;
S100:
    ix1 = ix+1;
    for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
    if(v <= f) goto S170;
    goto S30;
S130:
/*
     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
*/
    amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
    ynorm = -(k*k/(2.0*xnpq));
    alv = log(v);
    if(alv < ynorm-amaxp) goto S170;
    if(alv > ynorm+amaxp) goto S30;
/*
     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
     THE FINAL ACCEPTANCE/REJECTION TEST
*/
    x1 = ix+1.0;
    f1 = fm+1.0;
    z = n+1.0-fm;
    w = n-ix+1.0;
    z2 = z*z;
    x2 = x1*x1;
    f2 = f1*f1;
    w2 = w*w;
    if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
      (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
    goto S30;
S140:
/*
     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
*/
    qn = pow(q,(double)n);
    r = p/q;
    g = r*(n+1);
S150:
    ix = 0;
    f = qn;
    u = safe_ranf();
S160:
    if(u < f) goto S170;
    if(ix > 110) goto S150;
    u -= f;
    ix += 1;
    f *= (g/ix-r);
    goto S160;
S170:
    if(psave > 0.5) ix = n-ix;
    return ix;
}
int ignnbn(const int n,const double p)
/*
**********************************************************************
 
     int ignnbn(int n,double p)
                GENerate Negative BiNomial random deviate
                              Function
     Generates a single random deviate from a negative binomial
     distribution.
                              Arguments
     N  --> The number of trials in the negative binomial distribution
            from which a random deviate is to be generated.
	    JJV (N > 0)
     P  --> The probability of an event.
     JJV    (0.0 < P < 1.0)
                              Method
     Algorithm from page 480 of
 
     Devroye, Luc
 
     Non-Uniform Random Variate Generation.  Springer-Verlag,
     New York, 1986.
**********************************************************************
*/
{
static double y,a,r;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check Arguments
*/
    if(n <= 0L) ftnstop("N <= 0 in IGNNBN");
    if(p <= 0.0F) ftnstop("P <= 0.0 in IGNNBN");
    if(p >= 1.0F) ftnstop("P >= 1.0 in IGNNBN");
/*
     Generate Y, a random gamma (n,(1-p)/p) variable
     JJV Note: the above parametrization is consistent with Devroye,
     JJV       but gamma (p/(1-p),n) is the equivalent in our code
*/
    r = (double)n;
    a = p/(1.0F-p);
/*
 * JJV changed this to call SGAMMA directly
 *  y = gengam(a,r); <- OLD
 */
    y = sgamma(r)/a;
/*
     Generate a random Poisson(y) variable
*/
    return ignpoi(y);
}
int ignpoi(const double mu)
/*
**********************************************************************
     int ignpoi(double mu)
                    GENerate POIsson random deviate
                              Function
     Generates a single random deviate from a Poisson
     distribution with mean MU.
                              Arguments
     mu --> The mean of the Poisson distribution from which
            a random deviate is to be generated.
	    (mu >= 0.0)
     ignpoi <-- The random deviate.
                              Method
     Renames KPOIS from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Generation of Poisson Deviates
               From Modified Normal Distributions.
               ACM Trans. Math. Software, 8, 2
               (June 1982),163-179
**********************************************************************
**********************************************************************
                                                                      
                                                                      
     P O I S S O N  DISTRIBUTION                                      
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER GENERATION OF POISSON DEVIATES                
               FROM MODIFIED NORMAL DISTRIBUTIONS.                    
               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 
                                                                      
     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  
                                                                      
**********************************************************************
      INTEGER FUNCTION IGNPOI(IR,MU)
     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
             MU=MEAN MU OF THE POISSON DISTRIBUTION
     OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
     SEPARATION OF CASES A AND B
*/
{
static double a0 = -0.5;
static double a1 = 0.3333333;
static double a2 = -0.2500068;
static double a3 = 0.2000118;
static double a4 = -0.1661269;
static double a5 = 0.1421878;
static double a6 = -0.1384794;
static double a7 = 0.125006;
/* JJV changed the initial values of MUPREV and MUOLD */
static double muold = -1.0E37;
static double muprev = -1.0E37;
static double fact[10] = {
    1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
};
/* JJV added ll to the list, for Case A */
 static int ipoi,j,k,kflag,l,ll,m;
static double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,
    t,u,v,x,xx,pp[35];

    if(mu == muprev) goto S10;
    if(mu < 10.0) goto S120;
/*
     C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
     JJV changed l in Case A to ll
*/
    muprev = mu;
    s = sqrt(mu);
    d = 6.0*mu*mu;
/*
             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
             PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
*/
    ll = (int) (mu-1.1484);
S10:
/*
     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
*/
    g = mu+s*snorm();
    if(g < 0.0) goto S20;
    ipoi = (int) (g);
/*
     STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
*/
    if(ipoi >= ll) return ipoi;
/*
     STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
*/
    fk = (double)ipoi;
    difmuk = mu-fk;
    u = safe_ranf();
    if(d*u >= difmuk*difmuk*difmuk) return ipoi;
S20:
/*
     STEP P. PREPARATIONS FOR STEPS Q AND H.
             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
*/
    if(mu == muold) goto S30;
    muold = mu;
    omega = 0.3989423/s;
    b1 = 4.166667E-2/mu;
    b2 = 0.3*b1*b1;
    c3 = 0.1428571*b1*b2;
    c2 = b2-15.0*c3;
    c1 = b1-6.0*b2+45.0*c3;
    c0 = 1.0-b1+3.0*b2-15.0*c3;
    c = 0.1069/mu;
S30:
    if(g < 0.0) goto S50;
/*
             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
*/
    kflag = 0;
    goto S70;
S40:
/*
     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
*/
    if(fy-u*fy <= py*exp(px-fx)) return ipoi;
S50:
/*
     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
*/
    e = sexpo();
    u = safe_ranf();
    u += (u-1.0);
    t = 1.8+fsign(e,u);
    if(t <= -0.6744) goto S50;
    ipoi = (int) (mu+s*t);
    fk = (double)ipoi;
    difmuk = mu-fk;
/*
             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
*/
    kflag = 1;
    goto S70;
S60:
/*
     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
*/
    if(c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) goto S50;
    return ipoi;
S70:
/*
     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
             CASE IPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
*/
    if(ipoi >= 10) goto S80;
    px = -mu;
    py = pow(mu,(double)ipoi)/ *(fact+ipoi);
    goto S110;
S80:
/*
             CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
             A0-A7 FOR ACCURACY WHEN ADVISABLE
             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
*/
    del = 8.333333E-2/fk;
    del -= (4.8*del*del*del);
    v = difmuk/fk;
    if(fabs(v) <= 0.25) goto S90;
    px = fk*log(1.0+v)-difmuk-del;
    goto S100;
S90:
    px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
S100:
    py = 0.3989423/sqrt(fk);
S110:
    x = (0.5-difmuk)/s;
    xx = x*x;
    fx = -0.5*xx;
    fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
    if(kflag <= 0) goto S40;
    goto S60;
S120:
/*
     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
     JJV changed MUPREV assignment to initial value
*/
    muprev = -1.0E37;
    if(mu == muold) goto S130;
/* JJV added argument checker here */
    if(mu >= 0.0) goto S125;
    (void)fprintf(stderr,"MU < 0 in IGNPOI: MU %16.6E\n",mu);
    (void)fputs("Abort\n",stderr);
    exit(EXIT_FAILURE);
S125:
    muold = mu;
    m = max(1,(int) (mu));
    l = 0;
    p = exp(-mu);
    q = p0 = p;
S130:
/*
     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
*/
    u = safe_ranf();
    ipoi = 0;
    if(u <= p0) return ipoi;
/*
     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
             (0.458=PP(9) FOR MU=10)
*/
    if(l == 0) goto S150;
    j = 1;
    if(u > 0.458) j = min(l,m);
    for(k=j; k<=l; k++) {
        if(u <= *(pp+k-1)) goto S180;
    }
    if(l == 35) goto S130;
S150:
/*
     STEP C. CREATION OF NEW POISSON PROBABILITIES P
             AND THEIR CUMULATIVES Q=PP(K)
*/
    l += 1;
    for(k=l; k<=35; k++) {
        p = p*mu/(double)k;
        q += p;
        *(pp+k-1) = q;
        if(u <= q) goto S170;
    }
    l = 35;
    goto S130;
S170:
    l = k;
S180:
    return k;
}
int ignuin(const int low,const int high)
/*
**********************************************************************
     int ignuin(int low,int high)
               GeNerate Uniform INteger
                              Function
     Generates an integer uniformly distributed between LOW and HIGH.
                              Arguments
     low --> Low bound (inclusive) on integer value to be generated
     high --> High bound (inclusive) on integer value to be generated
                              Note
     If (HIGH-LOW) > 2,147,483,561 prints error message on * unit and
     stops the program.
**********************************************************************
     IGNLGI generates integers between 1 and 2147483562
     MAXNUM is 1 less than maximum generable value
*/
{
static int ranp1;
	
    if(!(low > high)) goto S10;
    (void)fputs(" low > high in ignuin - ABORT\n",stderr);
    exit(EXIT_FAILURE);

S10:

    if(!(low == high)) goto S30;
    return low;

S30:
    return low+(int)(safe_ranf()*(double)(ranp1+1));
}

int mltmod(const int a,const int s,const int m)
/*
**********************************************************************
     int mltmod(int a,int s,int m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
{
#define h2 15
#define h (1<<h2)
	
static int a0,a1,k,p,q,qh,rh;
/*
     H2 = ((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
    if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
    (void)fputs(" a, m, s out of order in mltmod - ABORT!\n",stderr);
    (void)fprintf(stderr," a = %12d s = %12d m = %12d\n",a,s,m);
    (void)fputs(" mltmod requires: 0 < a < m; 0 < s < m\n",stderr);
    exit(EXIT_FAILURE);
S10:
    if(!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
S20:
	a1 = a>>h2;
	  a0 = a-(a1<<h2);
    qh = m>>h2;
	  rh = m-(qh<<h2);
    if(!(a1 >= h)) goto S50;
    a1 -= h;
    k = s/qh;
    p = (int)(h*(s-k*qh)-k*rh);
S30:
    if(!(p < 0)) goto S40;
    p += m;
    goto S30;
S40:
    goto S60;
S50:
    p = 0;
S60:
/*
     P = (A2*S*H)MOD M
*/
    if(!(a1 != 0)) goto S90;
    q = m/a1;
    k = s/q;
    p -= (k*(m-a1*q));
    if(p > 0) p -= m;
    p += (a1*(s-k*q));
S70:
    if(!(p < 0)) goto S80;
    p += m;
    goto S70;
S90:
S80:
    k = p/qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
    p = (int)(h*(p-k*qh)-k*rh);
S100:
    if(!(p < 0)) goto S110;
    p += m;
    goto S100;
S120:
S110:
    if(!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
    q = m/a0;
    k = s/q;
    p -= (k*(m-a0*q));
    if(p > 0) p -= m;
    p += (a0*(s-k*q));
S130:
    if(!(p < 0)) goto S140;
    p += m;
    goto S130;
S150:
S140:
    return p;
#undef h
}

double sexpo(void) /*@*/
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER METHODS FOR SAMPLING FROM THE                 
               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
{
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
    .9999999
};
static int i;
static double a,u,ustar,umin;
static double *q1 = q;
S21:
    a = 0.0;
    u = safe_ranf();
    goto S30;
S20:
    a += *q1;
S30:
    u += u;
/*
 * JJV changed the following to reflect the true algorithm and prevent
 * JJV unpredictable behavior if U is initially 0.5.
 *  if(u <= 1.0) goto S20;
 */
    if(u < 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    return a+u;
S60:
	 if(u>q[7]) goto S21;
    i = 1;
    ustar = safe_ranf();
    umin = ustar;
S70:
    ustar = safe_ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    return a+umin**q1;
}

double sgamma(const double a)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  G A M M A  DISTRIBUTION                             
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
               PARAMETER  A >= 1.0  !                                 
                                                                      
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               GENERATING GAMMA VARIATES BY A                         
               MODIFIED REJECTION TECHNIQUE.                          
               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  
                                                                      
     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
                                 (STRAIGHTFORWARD IMPLEMENTATION)     
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
                                                                      
               PARAMETER  0.0 < A < 1.0  !                            
                                                                      
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
               COMPUTING, 12 (1974), 223 - 246.                       
                                                                      
     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    
                                                                      
**********************************************************************
     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
*/
{
static double q1 = 4.166669E-2;
static double q2 = 2.083148E-2;
static double q3 = 8.01191E-3;
static double q4 = 1.44121E-3;
static double q5 = -7.388E-5;
static double q6 = 2.4511E-4;
static double q7 = 2.424E-4;
static double a1 = 0.3333333;
static double a2 = -0.250003;
static double a3 = 0.2000062;
static double a4 = -0.1662921;
static double a5 = 0.1423657;
static double a6 = -0.1367177;
static double a7 = 0.1233795;
static double e1 = 1.0;
static double e2 = 0.4999897;
static double e3 = 0.166829;
static double e4 = 4.07753E-2;
static double e5 = 1.0293E-2;
static double aa = 0.0;
static double aaa = 0.0;
static double sqrt32 = 5.656854;
/* JJV added b0 to fix rare and subtle bug */
 static double sgam,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p;
    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = snorm();
    x = s+0.5*t;
    sgam = x*x;
    if(t >= 0.0) return sgam;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = safe_ranf();
    if(d*u <= t*t*t) return sgam;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463+s+0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-1*s;
S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if(x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if(log(1.0-u) <= q) return sgam;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = sexpo();
    u = safe_ranf();
    u += (u-1.0);
    t = b+fsign(si*e,u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if(t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
/*
 * JJV modified the code through line 115 to handle large Q case
 */
    if(q < 15.0) goto S95;
/*
 * JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
 * JJV so reformulate test at 110 in terms of one EXP, if not too big
 * JJV 87.49823 is close to the largest real which can be
 * JJV exponentiated (87.49823 = log(1.0E38))
 */
    if((q+e-0.5*t*t) > 87.49823) goto S115;
    if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
    goto S115;
S95:
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
S115:
    x = s+0.5*t;
    return x*x;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))

     JJV changed B to B0 (which was added to declarations for this)
     JJV in 120 to END to fix rare and subtle bug.
     JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
     JJV Reasons: the state of AA only serves to tell the A >= 1.0
     JJV case if certain A-dependent constants need to be recalculated.
     JJV The A < 1.0 case (here) no inter changes any of these, and
     JJV the recalculation of B (which used to change with an
     JJV A < 1.0 call) is governed by the state of AAA anyway.
    aa = 0.0;
*/
    b0 = 1.0+0.3678794*a;
S130:
    p = b0*safe_ranf();
    if(p >= 1.0) goto S140;
    sgam = exp(log(p)/ a);
    if(sexpo() < sgam) goto S130;
    return sgam;
S140:
    sgam = -log((b0-p)/ a);
    if(sexpo() < (1.0-a)*log(sgam)) goto S130;
    return sgam;
}
double snorm(void)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  N O R M A L  DISTRIBUTION                           
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
               SAMPLING FROM THE NORMAL DISTRIBUTION.                 
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
{
static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static int i;
static double sn,u,s,ustar,aa,w,y,tt;
    u = safe_ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (int) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar > *(t+i-1)) {
		 w = (ustar-*(t+i-1))**(h+i-1);
		 goto S50;
	 }
/*
                                CENTER CONTINUED
*/
    u = safe_ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = safe_ranf();
S80:
    if(ustar > tt) goto S50;
    u = safe_ranf();
    if(ustar >= u) goto S70;
    ustar = safe_ranf();
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = safe_ranf();
    if(ustar <= tt) {
		 u = safe_ranf();
		 if(ustar >= u) goto S150;
		 u = safe_ranf();
		 goto S140;
	 }
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    sn = y;
    if(s == 1.0) sn = -y;
    return sn;
}
