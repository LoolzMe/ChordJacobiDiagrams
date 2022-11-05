	/* mod P arithmetic */
#ifndef _MODP_h
#define _MODP_h 1

#define P	97	// Guess what this is.
#include <stdio.h>

class MODP
{
private:
  short int x;
  short int normalize(short int i)
	{short int a; a=i%P; if (a<0) a=P+a; return a;}

public:

  MODP() {x=0;}
  MODP(short int i) :x(normalize(i)) {}

  //void operator = (const MODP& r) {x=r.x;}
  void operator += (const MODP& r) {x=normalize(x+r.x);}
  void operator *= (const MODP& r) {x=normalize(x*r.x);}

  int iszero() {return 0 == x;}
  int isnotzero() {return 0 != x;}

  int operator == (const MODP& r) {return x == r.x;}
  int operator != (const MODP& r) {return x != r.x;}
  int operator < (const MODP& r) {return x < r.x;}

  MODP operator + (const MODP& r) {return (MODP) (x+r.x);}
  MODP operator - () {return (MODP) (P-x);}
  MODP operator * (const MODP& r) {return (MODP) (x*r.x);}
  MODP operator / (const MODP& r)
  {
    short int quo,temp;
    short int prev=0; short int last=1; short int big=P, small=r.x;
    while (small!=1) {
      temp=small; small=big-small*(quo=big/small); big=temp;
      temp=prev; prev=last; last=temp-quo*last;
    }
    return (MODP) (x*last);
  }

  short int operator ! () {return (x>P/2 ? x-P : x);}
};

#endif
#define CPP
// MODP.h"
#define UNUSED(name)		name	/* An ugly hack */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define M 			5	/* No. of loops */
#define framed			0	/* framed flag: */
					/* 0 - consider only legal diagrams */
					/* 1 - consider all diagrams */
#define solveQ			1	/* solve? */
					/* 0 - no. */
					/* 1 - yes, consider both parities */
					/* 2 - yes, even parity */
					/* 3 - yes, odd parity */
#define statistics_interval	10	/* Statistics printing interval */
#define treat_vans		0	/* Simplify when a van is found? */
#define treat_eqls		0	/* Simplify when an eql is found? */
#define afs_depth		4	/* Depth of afs array */
#define afs_step		4	/* afs step size in bits */
#define cleanup_ratio		15	/* nfree shrinks to #%: do cleanup */

#define CASE(a1,a2,a3,a4,a5,a6,a7,a8,a9) \
					(M==1 ? a1 : (M==2 ? a2 : \
					(M==3 ? a3 : (M==4 ? a4 : \
					(M==5 ? a5 : (M==6 ? a6 : \
					(M==7 ? a7 : (M==8 ? a8 : \
					(M==9 ? a9 : 0 )))))))))
#if framed
#define ndiags	CASE(1,2,5,18,105,902,9749,0,0)
#else
#define ndiags	CASE(0,1,2,7,36,300,3218,42335,644808)
#endif

#define N			8	/* Maximal number of WSs to */
					/* participate in a single product */
#define max_C_bank		32	/* Maximal size of a C bank */
#define max_m			16	/* Maximal m is this program */
/* #define abs(n)		(n>0 ? n : -n)	The absolute value function */
#define max(a,b)	((a)>(b) ? (a) : (b))	/* Maximum */
#define min(a,b)	((a)<(b) ? (a) : (b))	/* Minimum */

		/* Data types */
#ifndef CPP
typedef int COEFF;
#else
typedef MODP COEFF;
#endif
struct ACCUMULATOR {
  int len;
  COEFF *coeff;
};
struct EQUATION {
  int len;
  int *at;
  COEFF *coeff;
};

		/* Procedures */
void generate(int);			/* Generate diagrams */
int solve();				/* Solve system */
void print_solution();			/* Final printing of a solution */
void print_diagram(int m,int count);	/* Print a diagram */
int cyclic_order();			/* Check if i,j,k are cyclically */
					/* ordered */
void make_v(int m,int ae[],int v[]);	/* Make v from ae */
int find();				/* Search in a lexicographicaly */
					/* oredered array */
void cross(int m,int ae[],int v[],int d,int C[]);
					/* The m==2 weight system */
void m_is_3(int m,int ae[],int v[],int d,int C[]);
					/* The m==3 weight system */
void su(int m,int ae[],int v[],int d,int C[]);
					/* Computation of the su(N) weights */
void su_hat(int m,int ae[],int v[],int d,int C[]);
					/* Computation of su_hat(N) weights */
void su_general(int m,int ae[],int v[],int d,int C[],int hatQ);
					/* Computation of the general (hat */
					/* and n-hat) su(N) weights */
void so(int m,int ae[],int v[],int d,int C[]);
					/* Computation of the so(N) weights */
void so_hat(int m,int ae[],int v[],int d,int C[]);
					/* Computation of so_hat(N) weights */
void sp(int m,int ae[],int v[],int d,int C[]);
					/* Computation of sp(N) weights */
void sp_hat(int m,int ae[],int v[],int d,int C[]);
					/* Computation of sp_hat(N) weights */
void product(int m,int ae[],int v[],int n,int m_[],
	void (*(ws_[]))(int m,int ae[],int v[],int d,int C[]),
	int d_[],int C[]);		/* The WS product computation */
					/* routine */
void adjoint(int m,int ae[],int v[],int d,
	void (*ws)(int m,int ae[],int v[],int d,int C[]),int C[]);
					/* Computations in the Adjoint */
					/* representation */
int hat();				/* `hatting' a representation */
void best_version(int m,int ae[]);	/* Find the best version of ae */
int best_versionQ(int m,int ae[]);	/* Is this the best version of ae? */
void init_solve();
void generate_ee();
void print_status();
void generate_equation();
void addto(struct ACCUMULATOR *acc, struct EQUATION eq);
                                        /* Add eq2 into eq1 with the */
                                        /* correct coefficient */
void new_entry();
void interchange();
struct EQUATION compress(struct ACCUMULATOR acc);
void new_equality(int at);
void new_van(int at);
void print_statistics();
int find_diag(unsigned int aec);	/* Find a diagram in diags */
void to_basis(int ebound,int ibound);	/* Express acc in terms of basis */
					/* but never let len>bound */
void printCD(int v[2*M]);		/* print in CD format */
void print_d(int m,int count);		/* Print in d format */
int next_ae(int m,int ab[],int ae[]);	/* Find next ae */
unsigned int compress_ae(int m,int ae[]);
					/* Compress ae */
void uncompress_ae(int m,unsigned int aec,int ae[]);
                                        /* Uncompress ae */
void do_parity();			/* fill matrix with either odd or */
					/* even layout */
void afs_update(int position, int dal);	/* update afs by dal at position */
int afs_search(int from);		/* Search for the next non-zero */
					/* entry in acc, past&incl. from. */
void cleanup();				/* Upward simplify matrix */
int statisticsQ();			/* Time to print statistics? */

		/* Large arrays */
unsigned int diags[ndiags];		/* ae of the diagrams */
struct EQUATION matrix[ndiags];		/* The Vassiliev equations */

		/* Smaller arrays */
int binomial[M+1][M+1];			/* The binomial coefficients */
extern int M2;				/* No. of vertices */
extern int count;			/* No. of diagrams */
int rep;				/* Current report number */

		/* solve.c variables */
struct ACCUMULATOR acc;			/* The equation currently in use */
int *(afs[afs_depth]);				/* acc fill summary */
int ae[2*M];				/* The other end of the arc */
					/* connected to the current vertex */
int v[2*M];				/* The vertices */
int ee[2*M-2];				/* The equation symbol */
int i,j,k,l,q;				/* Temporary variables */
int ca;					/* Current arc */
int cab,cae;				/* Current arc beginning and end */
int sa,side,sab,sae,sa_ends[2];		/* Special arc in eq. generator */
int tof,tofree[M];			/* Which free vertex to connect next */
int after,shift;			/* Used in the eq. construction */
int pivot;				/* The pivot in the current eqn.  */
int max_len,max_length;			/* Length of the current eqn.  */
int s,t;				/* Euclid multipliers */
int deg;				/* No. of degrees of freedom */
static int nfree=ndiags;		/* No. of free rows in matrix */
int eqs,ads;				/* Statistics */
int iteration,no_change,last_new,
      all_done;				/* iteration control */
int intrchngs,failures;			/* Insertion statistics */
int tl;					/* Total length of `matrix' */
int from,to;				/* Searching domain */
int position;				/* Position to which a new diagram */
					/* should be inserted */
struct EQUATION term;			/* The new term while equation */
					/* is being generated. */
int eqls;				/* Number of equalities encountered */
int vans;				/* Number of zeros encountered */
int next_cleanup;			/* time (in nfree) of next cleanup */
// ws.h"

void init();			/* Initializations */

int M2=2*M;
int count;

main()
{
  int all_done=0;		/* Whether solve() has succeded */
  init();
  generate(M);
  if (solveQ) {
    all_done=solve();
    if (!all_done) printf("\n\n(*	WARNING:  all_done flag is off!!!!	*)\n\n");
  }
  print_solution();
}

void init()
{

  int i,j;				 /* Temporary variables */
  int ca=0;				/* Current arc */

		/* Initializations */
  printf("(*	The degree %d case: solveQ=%d, P=%d. *)\n",M,solveQ,P);
  rep=0;
  for (i=0; i<=M; ++i) {
    binomial[i][0]=1; binomial[i][i]=1;
    for (j=1; j<i; ++j) binomial[i][j]=binomial[i-1][j-1]+binomial[i-1][j];
  }
}
// ws.h"
int arc[2*M];			/* Current arc number */

void generate(int m)
{
  int ae[2*M];			/* The other end of the arc */
				/* connected to the current vertex */
  int ab[M];			/* Beginning of the current arc */
  int count;			/* # of current diagram */
  unsigned int aec;		/* ae compressed */
  int i;			/* Temporary variable */

		/* Initializing count,ab,ae */
  count=0;
  if (framed) for (i=0; i<m; ++i) {
    arc[2*i]=i; arc[2*i+1]=i;
    ab[i]=2*i; ae[2*i]=2*i+1; ae[2*i+1]=2*i;
  } else {
    for (i=0; i<m/2; ++i) {
      ab[2*i]=4*i; ab[2*i+1]=4*i+1;
      arc[4*i]=2*i; arc[4*i+2]=2*i; ae[4*i]=4*i+2; ae[4*i+2]=4*i;
      arc[4*i+1]=2*i+1; arc[4*i+3]=2*i+1; ae[4*i+1]=4*i+3; ae[4*i+3]=4*i+1;
    }
    if (m%2!=0) {
      arc[2*m-2]=m-2; ae[2*m-5]=2*m-2; ae[2*m-2]=2*m-5;
      arc[2*m-3]=m-1; arc[2*m-1]=m-1;
      ab[m-1]=2*m-3; ae[2*m-3]=2*m-1; ae[2*m-1]=2*m-3;
    }
  }

		/* Main loop */
  do {
    if (statisticsQ()) {printf("generating %d*)\n",count); fflush(stdout);}
    diags[count++]=compress_ae(m,ae);
  } while (next_ae(m,ab,ae));

		/* Reversing the order */
  --count;
  for (i=0; i<=count/2; ++i)
  {
    if (statisticsQ()) {printf("reversing %d*)\n",i); fflush(stdout);}
    aec=diags[i]; diags[i]=diags[count-i]; diags[count-i]=aec;
  }
  ++count;

		/* Printing */
  if (framed) printf("\nDimCD[%d]=%d\n\n",M,count);
  else printf("\nDimCDr[%d]=%d\n\n",M,count);
}

int next_ae(int m,int ab[],int ae[])
{
  int m2;			/* No. of vertices */
  int at;			/* try pushing this arc */
  int cab;			/* Current arc beginning */
  int foundQ;			/* good new ae found? */

  m2=2*m;
  at=m-2;
  while (at<m && at>=0) {
    foundQ=0;
    for (i=1+ae[cab=ab[at]]; !foundQ && i<m2; ++i) {
      if (arc[i]>=at) foundQ=1;
      if (foundQ && i-cab < ae[0]) foundQ=0;
      if (foundQ && !framed && 1==(i-cab)) foundQ=0;
      if (foundQ && cab!=ae[cab]) arc[ae[cab]]=m;
      if (foundQ) {
        if (cab!=ae[cab]) arc[ae[cab]]=m;
        arc[i]=at; ae[cab]=i; ae[i]=cab;
        if (at==m-1 && !best_versionQ(m,ae)) foundQ=0;
      }
    } 
    at+=foundQ;
    if (foundQ && m>at) {
      i=cab+1;
      while (at>arc[i]) ++i;
      ab[at]=i; arc[i]=at; ae[i]=i;
    } else if (!foundQ) {
      arc[cab]=m; arc[ae[cab]]=m;
      if (0==(--at) && ae[0]>m) --at;
    }
  }
  if (at==m) return(1); else return(0);
}

unsigned int compress_ae(int m,int ae[])
{
  unsigned int aec;			/* ae compressed */
  int i,j;				/* Temporary */
  int m2;				/* No. of vertices */

  m2=2*m; aec=0; j=0;
  for (i=0; i<m2; ++i) if (j<m-1 && ae[i]>i) {
    aec=16*aec+(ae[i]-i-1);
    ++j;
  }
  return(aec);
}

void uncompress_ae(int m,unsigned int aec,int ae[])
{
  unsigned int taec;                    /* temporary ae compressed */
  int diff[M];				/* aec in base 16 */
  int i,j;                              /* Temporary */
  int m2;                               /* No. of vertices */

  m2=2*m; taec=aec;
  for (i=m-2; i>=0; --i) {
    diff[i]=taec%16+1;
    taec/=16;
  }
  for (i=0; i<m2; ++i) ae[i]=m2;
  j=0;
  for (i=0; i<m-1; ++i) {
    while (ae[j]<m2) ++j;
    ae[j]=j+diff[i];
    ae[ae[j]]=j;
  }
  while (ae[j]<m2) ++j;
  i=j+1;
  while (ae[i]<m2) ++i;
  ae[j]=i; ae[i]=j;
}
// ws.h"

int solve()
{
  int l1,l2,l3;		/* Lengths of segments in ee */
  init_solve();

  while (!all_done && solveQ) {
    if (solveQ>1) do_parity();

    next_cleanup=(cleanup_ratio*nfree)/100;
    if (next_cleanup>0) printf("(*First cleanup will be at nfree=%d.*)\n",
								next_cleanup);
    eqs=last_new=0;
    all_done=no_change=1;
    for (i=0; i<M-1; ++i) tofree[i]=0;

		/* The equations loop */
    for (sa=0; sa<M-1; ++sa) do {
      generate_ee();
      if (!framed) for (i=0; i<M2-2; ++i) if (ee[i]-i == 1) goto irrelevant;
      if (solveQ>1) {
        l1=sa_ends[1]-sa_ends[0];
        l2=M2-2-sa_ends[1];
        l3=1+sa_ends[0];
        if (((l1>l2 || l2>l3) && (l1<l2 || l2<l3))) goto irrelevant;
      }
      if (statisticsQ()) print_status();
      if (nfree<=next_cleanup && nfree>0) cleanup();
      ++eqs;
      generate_equation();
      while (acc.len>0) {
        if (acc.len>max_len) max_len=acc.len;
        pivot=afs_search(pivot);
        if (!matrix[pivot].len) new_entry();
        else {
          interchange();
          ++ads;
          addto(&acc,matrix[pivot]);
        }
      }
irrelevant:
      i=M-3; j=0;
      while (tofree[i]==(j+=2) && i>=0) tofree[i--]=0;
      if (i>=0) ++tofree[i];
    } while (nfree>0 && i>=0);
  }
  return(all_done);
}

void addto(struct ACCUMULATOR *acc, struct EQUATION eq)
{
  COEFF coeff;				/* Temporary variables */
  COEFF factor;				/* Numerator and denumerator of */
					/* factor currently in use. */
  int i;                             	/* Temporary variable */
  int dal;				/* delta acc length */

  if (!eq.len) return;
  factor= -acc->coeff[eq.at[0]]/eq.coeff[0];
  for (i=0; i<eq.len; ++i) {
    coeff=acc->coeff[eq.at[i]];
    if (coeff.isnotzero()) dal=-1; else dal=0;
    coeff+=factor*eq.coeff[i];
    if (coeff.isnotzero()) ++dal;
    acc->len+=dal; if (dal) afs_update(eq.at[i],dal);
    acc->coeff[eq.at[i]]=coeff;
  }
}

void init_solve()
{
  int i,j,k;				/* Temporary variables */

  max_len=0 ; max_length=0; eqs=ads=0;
  no_change=last_new=all_done=0;
  intrchngs=failures=0; eqls=0; vans=0;
  nfree=ndiags;
  for (count=0; count<ndiags; ++count) matrix[count].len=0;
  acc.coeff=new (COEFF [ndiags]);
  acc.len=0; for(i=0; i<ndiags; ++i) acc.coeff[i]=0;
  k=ndiags;
  for (i=0; i<afs_depth; ++i) {
    k>>=afs_step;
    afs[i]=new (int [k+1]);
    for (j=0; j<=k; ++j) afs[i][j]=0;
  }
}

void generate_ee()
{
  for (i=0; i<M2-2; ++i) ee[i]= -1;
  for (ca=0; ca<M-1; ++ca) {
    tof=tofree[ca];
    for (cab=ca; ee[cab]>=0; ++cab);
    cae=cab;
    for (i=0; i<=tof; ++i) {
      ++cae;
      while (ee[cae]>=0) ++cae;
    }
    ee[cab]=cae; ee[cae]=cab;
    if (ca==sa) {sa_ends[0]=cab; sa_ends[1]=cae;}
  }
}

void print_status()
{
  int count;				/* Temporary variables */

  tl=0;
  for (count=0; count<ndiags; ++count) tl+=matrix[count].len;
  //printf("(*eqs=%d tl=%d ads=%d intrchngs=%d eqls=%d vans=%d nfree=%d*)\n",
	  //eqs,tl,ads,intrchngs,eqls,vans,nfree);
  printf("eqs=%d tl=%d ads=%d intrchngs=%d nfree=%d*)\n",
		eqs,tl,ads,intrchngs,nfree);
  fflush(stdout);
}

		/* Generating the equation */
void generate_equation()
{
  COEFF coeff;
  int dal;				/* delta acc length */
  unsigned int aec;			/* ae compressed */

  pivot=ndiags;
  for (side=0; side<2; ++side) for (after=0; after<2; ++after) {
    sab=sa_ends[side]; sae=sa_ends[1-side]; shift=0;
    for (i=0; i<M2-2; ++i) if (i==sab) {
      ae[M2-1]=sab+after;
      ae[sab+after]=M2-1;
      ae[sab+1-after]=ee[sab]+1-side;
      shift=1;
    } else ae[i+shift]=ee[i]+(ee[i]>sab);
    ae[sae+1-side]+=1-after;

    best_version(M,ae);
    aec=compress_ae(M,ae);
    position=find_diag(aec);

    if (position<ndiags && position>=0 && matrix[position].len!=1) {
      if (matrix[position].len==2) {
        coeff=-matrix[position].coeff[1]/matrix[position].coeff[0];
        position=matrix[position].at[1];
      } else coeff=1;
      pivot=min(pivot,position);
      if (acc.coeff[position].isnotzero()) dal=-1; else dal=0;
      acc.coeff[position]+=(after ? coeff : -coeff);
      if (acc.coeff[position].isnotzero()) ++dal;
      acc.len+=dal; if (dal) afs_update(position,dal);
    }
  }
}

void new_entry()
{
  int i;

  matrix[pivot]=compress(acc);
  if (0 != matrix[pivot].coeff) {
    acc.len=0;
    for (i=0; i<matrix[pivot].len; ++i) {
      acc.coeff[matrix[pivot].at[i]]=0;
      afs_update(matrix[pivot].at[i],-1);
    }
    last_new=eqs;
    --nfree;
    if (treat_vans && 1==matrix[pivot].len) new_van(pivot);
    if (treat_eqls && 2==matrix[pivot].len) new_equality(pivot);
  } else {
    printf("(*	OUT OF MEMORY!!	*)\n");
    all_done=0;
    print_statistics();
  }
}

void interchange()
{
  struct EQUATION mat;

  if (matrix[pivot].len<=acc.len) return;
  mat=matrix[pivot];
  matrix[pivot]=compress(acc);
  if (0 != matrix[pivot].coeff) {
    for (i=0; i<matrix[pivot].len; ++i) {
      acc.coeff[matrix[pivot].at[i]]=0;
      afs_update(matrix[pivot].at[i],-1);
    }
    acc.len=mat.len;
    for (i=0; i<acc.len; ++i) {
      acc.coeff[mat.at[i]]=mat.coeff[i];
      afs_update(mat.at[i],1);
    }
#ifndef CPP
    free((char *) mat.at);
    free((char *) mat.coeff);
#else
    delete mat.at;
    delete mat.coeff;
#endif
    no_change=0;
    ++intrchngs;
    if (treat_vans && 1==matrix[pivot].len) new_van(pivot);
    if (treat_eqls && 2==matrix[pivot].len) new_equality(pivot);
  } else {
    printf("(*	OUT OF MEMORY!!	*)\n");
    all_done=0;
    print_statistics();
  }
}

struct EQUATION compress(struct ACCUMULATOR acc)
{
  struct EQUATION eq;		/* Temporary storage for compressed acc */
  int at;			/* Current at */
  int i;			/* Temporary variables */

  eq.at=new (int [acc.len]);
  eq.coeff=new (COEFF [acc.len]);
  eq.len=acc.len;

  i=0; at=afs_search(0);
  while (at<ndiags) {
    eq.at[i]=at;
    eq.coeff[i]=acc.coeff[at];
    ++i; at=afs_search(at+1);
  }
  return(eq);
}

void new_equality(int at)
{
  int change_to,count;
  int i,saved;
  COEFF coeff,put;
  struct EQUATION mat;

  ++eqls;
  change_to=matrix[at].at[1];
  coeff=-matrix[at].coeff[1]/matrix[at].coeff[0];
  for (count=0; count<at; ++count) if (matrix[count].len>1) {
    saved=0;
    mat=matrix[count];
    for (i=1; mat.at[i]< at && i<mat.len; ++i);
    if (i<mat.len && mat.at[i]==at) {
      ++saved;
      put=coeff*mat.coeff[i];
      ++i;
      while (mat.at[i]<change_to && i<mat.len) {
        mat.at[i-1]=mat.at[i]; mat.coeff[i-1]=mat.coeff[i];
        ++i;
      }
      if (i<mat.len) if (mat.at[i]==change_to) {
        mat.at[i-1]=change_to;
        mat.coeff[i-1]=put+mat.coeff[i];
        if (mat.coeff[i-1].isnotzero()) {
          ++i;
          while (i<mat.len) { 
            mat.at[i-1]=mat.at[i]; mat.coeff[i-1]=mat.coeff[i];
            ++i; 
          }
        } else {
          ++saved; ++i; 
          while (i<mat.len) {  
            mat.at[i-2]=mat.at[i]; mat.coeff[i-2]=mat.coeff[i]; 
            ++i;  
          }
        }
      } else {
        --saved;
        mat.at[i-1]=change_to; mat.coeff[i-1]=put;
      } else {
        --saved; 
        mat.at[i-1]=change_to; mat.coeff[i-1]=put;
      }
      if (saved) {
        mat.len-=saved;
        matrix[count].at=new (int [mat.len]);
        matrix[count].coeff=new (COEFF [mat.len]);
        for (i=0; i<mat.len; ++i) {
          matrix[count].at[i]=mat.at[i];
          matrix[count].coeff[i]=mat.coeff[i];
        }
        matrix[count].len=mat.len;
        delete mat.at;
        delete mat.coeff;
        if (treat_vans && mat.len==1) new_van(count);
        if (treat_eqls && mat.len==2) new_equality(count);
      }
    }
  }
}

int best_versionQ(int m,int ae[])
{
  int i,q;				/* Temporary */
  int m2;                               /* No. of vertices */
  int ret;				/* return value */

  m2=2*m;
  ret=1;
  for (q=1; q<m2 && ret; ++q) {
    for (i=0; ae[i]==(m2+ae[(i+q)%m2]-q)%m2 && i<m2-1; ++i);
    if (ae[i]>(m2+ae[(i+q)%m2]-q)%m2) ret=0;
  }
  return(ret);
}

void best_version(int m,int ae[])
{
  int temp[2*M];		/* temporary `ae' */
  int i,q;			/* Temporary */
  int bq;			/* Best q so far */
  int m2;			/* No. of vertices */

  m2=2*m;
  bq=0;
  for (q=1; q<m2; ++q) {
    for (i=0; (m2+ae[(i+bq)%m2]-bq)%m2==(m2+ae[(i+q)%m2]-q)%m2 && i<m2-1; ++i);
    if ((m2+ae[(i+bq)%m2]-bq)%m2>(m2+ae[(i+q)%m2]-q)%m2) bq=q;
  }
  if (bq!=0) {
    for (i=0; i<m2; ++i) temp[i]=ae[i];
    for (i=0; i<m2; ++i) ae[i]=(m2+temp[(i+bq)%m2]-bq)%m2;
  }
}

void print_statistics()
{
  int i;

  if (!all_done)
      printf("\n\n(*	WARNING:  all_done flag is off!!!!	*)\n\n");
  for (i=0; i<ndiags; ++i) if (matrix[i].len>max_length)
      max_length= matrix[i].len;
  printf("(*	Total number of equations considered was:%6d.	*)\n",eqs);
  printf("(*	Total number of row operations performed was:%7d.	*)\n",ads);
  printf("(*	Maximal equation length is:%4d.	*)\n",max_length);
  printf("(*	Maximal equation length considered was:%4d.	*)\n",max_len);
  printf("(*	Total length of the equations stored:%8d.	*)\n",tl);
}

int find_diag(unsigned int aec)		/* Find a diagram in diags */
{
  int i;				/* Temporary variable */
  int from,to,position;			/* Search outcome */

  from=0; to=ndiags-1;
  while (to>from) {
    position=(from+to)/2;
    if (aec==diags[position]) from=to=position;
       else if (aec>diags[position]) to=--position;
               else from=++position;
  }
  if (position<ndiags && position>=0)
     if (aec!=diags[position]) position=ndiags;
  return(position);
}

void to_basis(int ebound, int ibound)
{
  int pivot;
  int behind;		/* Number of non-zero acc entries left behind */

  pivot=0; behind=0;
  while (behind<=ibound && acc.len>0 && acc.len<=ebound && pivot<ndiags) {
    pivot=afs_search(pivot);
    if (pivot<ndiags) {
      if (matrix[pivot].len) addto(&acc,matrix[pivot]);
        else {++behind; if (ndiags==++pivot) return;}
    }
  }
}

void new_van(int at)
{
  int count,i;			/* Temporaries */
  struct EQUATION mat;

  ++vans;
  for (count=0; count<at; ++count) if (matrix[count].len>1) {
    mat=matrix[count];
    for (i=1; mat.at[i]< at && i<mat.len; ++i);
    if (i<mat.len && mat.at[i]==at) {
      ++i;
      while (i<mat.len) {
        mat.at[i-1]=mat.at[i]; mat.coeff[i-1]=mat.coeff[i];
        ++i;
      }
      --mat.len;
      matrix[count].at=new (int [mat.len]);
      matrix[count].coeff=new (COEFF [mat.len]);
      for (i=0; i<mat.len; ++i) {
        matrix[count].at[i]=mat.at[i];
        matrix[count].coeff[i]=mat.coeff[i];
      }
      matrix[count].len=mat.len;
      delete mat.at;
      delete mat.coeff;
      if (treat_vans && mat.len==1) new_van(count);
      if (treat_eqls && mat.len==2) new_equality(count);
    }
  }
}

		/* Doing parity right */
void do_parity()
{
  int count,position;			/* positions in diags[] */
  int j,k;                              /* Temporary */
  int m2;                               /* No. of vertices */

  m2=2*M;
  for (count=0; count<ndiags; ++count) if (matrix[count].len==0) {
    if (statisticsQ()) {printf("casting %d*)\n",count); fflush(stdout);}
    uncompress_ae(M,diags[count],ae);
    for (j=0; j<M; ++j) {k=ae[j]; ae[j]=m2-1-ae[m2-1-j]; ae[m2-1-j]=m2-1-k;}
    best_version(M,ae);
    position=find_diag(compress_ae(M,ae));
    if (position==count && solveQ==3 && P>2) {
      matrix[count].len=1;
      matrix[count].at=new (int[1]);
      matrix[count].coeff=new (COEFF[1]);
      matrix[count].at[0]=count;
      matrix[count].coeff[0]=1;
      --nfree; ++vans;
    }
    if (position>count) {
      matrix[count].len=2; 
      matrix[count].at=new (int[2]);   
      matrix[count].coeff=new (COEFF[2]);
      matrix[count].at[0]=count;
      matrix[count].coeff[0]=1;
      matrix[count].at[1]=position;
      matrix[count].coeff[1]=(solveQ==2 ? -1 : 1);
      --nfree; ++eqls;
    }
  }
}

void afs_update(int position, int dal)
{
  int i,p;

  p=(position>>afs_step);
  for (i=0; i<afs_depth; ++i) {
    afs[i][p]+=dal;
    p>>=afs_step;
  }
}

int afs_search(int from)
{
  int ret;			/* return value */
  int l;			/* current level of search */
  int i;			/* current search position */

  ret=from;
  while (ret<ndiags && 0!=afs[0][ret>>afs_step] && acc.coeff[ret].iszero())
	++ret;
  if (ret==ndiags || acc.coeff[ret].isnotzero()) return(ret);
  else {
    l=0;
    ret>>=afs_step;
    while(l>=0 && ret<=(ndiags>>(afs_step*(l+1)))) {
      if (afs[l][ret]!=0) {--l; ret<<=afs_step;}
      else if (l<afs_depth-1 && afs[l+1][ret>>afs_step]==0)
	      {++l; ret>>=afs_step;}
      else ++ret;
    }
    if (l>=0) return(ndiags);
    else {
      while (ret<ndiags && acc.coeff[ret].iszero()) ++ret;
      return(ret);
    }
  }
}

void cleanup()
{
  int count;		/* Current equation cleaned */
  int improved;		/* Number of equations improved */
  struct EQUATION eq;	/* Temporary equation storage */
  int i;		/* Temporary */

  improved=0;
  printf("\n(*cleanup called at:*)\n(*");
  print_status();
  for (count=ndiags-2; count>=0; --count) {
    if (statisticsQ()) {printf("count=%d*)\n",count); fflush(stdout);}
    if (matrix[count].len>1) {
      acc.len=1; afs_update(count,1);
      acc.coeff[count]=-1;
      to_basis((matrix[count].len==2) ? 1 : ndiags, matrix[count].len);
      if (acc.len<matrix[count].len) {
        delete matrix[count].at;
        delete matrix[count].coeff;
        ++acc.len; afs_update(count,1); 
        acc.coeff[count]=1;
        matrix[count]=compress(acc);
        ++improved;
      } else {
        eq.at=new (int [matrix[count].len]);
        eq.coeff=new (COEFF [matrix[count].len]);
        eq.len=matrix[count].len;
        for (i=0; i<eq.len; ++i) {
          eq.at[i]=matrix[count].at[i];
          eq.coeff[i]=matrix[count].coeff[i];
        }
        delete matrix[count].at;
        delete matrix[count].coeff;
        matrix[count].at=eq.at;
        matrix[count].coeff=eq.coeff;
      }
    }
    i=count;
    while ((i=afs_search(i))<ndiags) {acc.coeff[i]=0; afs_update(i,-1); ++i;}
    acc.len=0;
  }
  next_cleanup=(cleanup_ratio*next_cleanup)/100;
  printf("(*%d equations improved, next cleanup at %d.*)\n\n",
						improved,next_cleanup);
  fflush(stdout);
}
// ws.h"

void make_v(int m,int ae[],int v[])
{
  int m2;		/* 2*m */
  int ca;		/* Current arc */
  int i;		/* Temporary variable */

  ca=0; m2=2*m;
  for (i=0; i<m2; ++i) if (i<ae[i]) v[i]=v[ae[i]]=(ca++);
}

	/* WS producing routines. */
	/* Format: WS(m,ae,v,d,C) computes d weight systems on the diagram */
	/*         presented by (m,ae,v) and stores the result on C. */

void cross(int UNUSED(m), int ae[],int UNUSED(v)[], int UNUSED(d), int C[])
{
  if (*ae==2) *C=1; else *C=0;
}

void m_is_3(int UNUSED(m), int ae[],int UNUSED(v)[], int UNUSED(d), int C[])
{
  if (ae[0]==3 && ae[1]==4 && ae[2]==5) *C=2;
  else if (ae[0]==3 && ae[1]==5 && ae[2]==4) *C=1;
  else if (ae[0]==4 && ae[1]==3 && ae[2]==5) *C=1;
  else if (ae[0]==2 && ae[1]==4 && ae[3]==5) *C=1;
  else *C=0;
}

		/* The group su(N) */
void su(int m,int ae[],int v[],int d,int C[])
{
  su_general(m,ae,v,d,C,0);
}

void su_hat(int m,int ae[],int v[],int d,int C[])
{
  su_general(m,ae,v,d,C,1);
}

void su_general(int m,int ae[],int v[],int d,int C[],int hatQ)
{
  int m2;			/* No. of vertices */
  int Csu[2*max_m+1];		/* su(N) weights */
  int open[max_m+1];		/* Binary counter */
  int visited[2*max_m+1];	/* Wilson loop segments visited */
  int seg;			/* Segment pointer. Points to the vertex */
				/* before the current segment. */
  int power,sign;		/* Power and sign of current cont. */
  int i;			/* Temp. var. */

  m2=2*m; d=min(d,(hatQ ? (m+3)/2 : 2*m+1));
  visited[m2]=1;
  for (i=0; i<d; ++i) Csu[i]=0;
  for (i=0; i<=m; ++i) open[i]=0;
  do {
    for (i=0; i<m2; ++i) visited[i]=0;
    power=(hatQ ? -1 : m-1);
    seg=0; sign=1;
    while (!visited[seg]) {
      ++power;
      while (!visited[seg]) {
        visited[seg]=1;
        if (open[v[seg=(seg+1) % m2]]) seg=ae[seg];
      }
      for (seg=1; seg<m2 && visited[seg]; ++seg) ;
    }
    for (seg=0; seg<m; ++seg) if (!open[seg]) {
      power+=(hatQ ? 1 : -1);
      sign= -sign;
    }
    Csu[power/2]+=sign;
    for (i=0; open[i]; ++i) open[i]=0;
    open[i]=1;
  } while (!open[m]);
  for (i=0; i<d; ++i) C[i]=Csu[i];
}

void so(int m,int ae[],int v[],int d,int C[])
{
  int m2;			/* No. of vertices */
  int Cso[(max_m+1)];		/* so(N) weights */
  int crossed[max_m+1];		/* Crossed bridge? */
  int visited[2*max_m+1];	/* Wilson loop segments visited */
  int seg;			/* segment pointer. Points to the vertex */
				/* before the current segment. */
  int power,sign;		/* Power and sign of current cont. */
  int forward;			/* Direction of travel */
  int i;			/* Temp. var. */

  m2=2*m;
  visited[m2]=1;
  for (i=0; i<=m; ++i) Cso[i]=0;
  for (i=0; i<=m; ++i) crossed[i]=0;
  do {
    for (i=0; i<m2; ++i) visited[i]=0;
    seg=power=0; forward=1; sign=1;
    while (!visited[seg]) {
      ++power;
      while (!visited[seg]) {
        visited[seg]=1;
        if (crossed[v[i=(seg+forward) % m2]]) {
          forward=!forward;
          seg=(m2+ae[i]-1+forward) % m2;
        } else seg=(m2+ae[i]-1+forward) % m2;
      }
      for (seg=1; seg<m2 && visited[seg]; ++seg) ;
    }

    for (seg=0; seg<m; ++seg) if (crossed[seg]) sign=-sign;
    Cso[--power]+=sign;
    for (i=0; crossed[i]; ++i) crossed[i]=0;
    crossed[i]=1;
  } while (!crossed[m]);
  for (i=0; i<min(d,m+1); ++i) C[i]=Cso[i];
}

void so_hat(int m,int ae[],int v[],int d,int C[])
{
  int m2;			/* No. of vertices */
  int Cso[(max_m+1)];		/* so_hat(N) weights */
  int arc_status[max_m+1];	/* Trinary counter */
				/* 0: regular bridge */
				/* 1: crossed bridge */
				/* 2: closed bridge */
  int visited[2*max_m+1];	/* Wilson loop segments visited */
  int seg;			/* segment pointer. Points to the vertex */
				/* before the current segment. */
  int power,sign,anomaly;	/* Power, sign, and anomaly of current cont. */
  int forward;			/* Direction of travel */
  int i;			/* Temp. var. */

  m2=2*m;
  visited[m2]=1;
  for (i=0; i<=m; ++i) Cso[i]=0;
  for (i=0; i<=m; ++i) arc_status[i]=0;
  do {
    for (i=0; i<m2; ++i) visited[i]=0;
    seg=power=anomaly=0; sign=forward=1;
    while (!visited[seg]) {
      ++power;
      while (!visited[seg]) {
        visited[seg]=1;
        switch (arc_status[v[i=(seg+forward) % m2]]) {
        case 0:
          seg=(m2+ae[i]-1+forward) % m2;
          break;
        case 1:
          forward=!forward;
          seg=(m2+ae[i]-1+forward) % m2;
          break;
        case 2:
          seg=(seg+(forward ? 1 : m2-1)) % m2;
          break;
        }
      }
      for (seg=1; seg<m2 && visited[seg]; ++seg) ;
    }

    for (seg=0; seg<m; ++seg) switch (arc_status[seg]) {
    case 0:
      break;
    case 1:
      sign= -sign;
      break;
    case 2:
      ++anomaly;
      break;
    }
    --power;
    for (i=0; i<=anomaly; ++i) {
      Cso[power+i]+=sign*binomial[anomaly][i];
      sign= -sign;
    }
    for (i=0; arc_status[i]==2; ++i) arc_status[i]=0;
    ++arc_status[i];
  } while (!arc_status[m]);
  for (i=0; i<min(d,m+1); ++i) C[i]=Cso[i];
}

void sp(int m,int ae[],int v[],int d,int C[])
{
  int m2;			/* No. of vertices */
  int Csp[max_m+1];		/* sp(N) weights */
  int invert[max_m+1];		/* Whether the current arc is PQ inverting */
  int type[max_m+1];		/* Surgery type on current arc */
  int isP[2*max_m];		/* Whether current segment is a P */
  int crossed[max_m];		/* Crossed arc? */
  int visited[2*max_m+1];	/* Wilson loop segments visited */
  int seg;			/* segment pointer. Points to the vertex */
				/* before the current segment. */
  int power,sign;		/* Power and sign of current cont. */
  int forward;			/* Direction of travel */
  int p;			/* P now? */
  int legit;			/* A legitimate P? */
  int i,j;			/* Temp. var. */

  m2=2*m;
  visited[m2]=1;
  for (i=0; i<=m; ++i) Csp[i]=0;
  for (i=0; i<=m; ++i) invert[i]=0;
  do {
    p=1;
    for (i=0; i<m2; ++i) isP[i]=p^=invert[v[i]];
    legit=1;
    for (i=1; i<m2; ++i) legit=legit && !(invert[v[i]] && isP[i]==isP[ae[i]]);
    if (legit) {
      for (i=0; i<=m; ++i) type[i]=0;
      do {
        visited[0]=seg=power=0; sign=1; forward=1;
        for (i=0; i<m2; ++i) {
          visited[i]=0;
          j=v[i];
          if (invert[j]) crossed[j]=type[j];
          else if (crossed[j]=(isP[i]!=isP[ae[i]])) if (isP[i]) sign=-sign;
        }
        while (!visited[seg]) {
          ++power;
          while (!visited[seg]) {
            visited[seg]=1;
            if (crossed[v[i=(seg+forward) % m2]]) {
              forward=!forward;
              seg=(m2+ae[i]-1+forward) % m2;
            } else seg=(m2+ae[i]-1+forward) % m2;
          }
          for (seg=1; seg<m2 && visited[seg]; ++seg) ;
        }

        Csp[(--power)]+=sign;
        for (i=0; i<m && (type[i] || !invert[i]); ++i) type[i]=0;
        type[i]=1;
      } while (!type[m]);
    }
    for (i=0; invert[i]; ++i) invert[i]=0;
    invert[i]=1;
  } while (!invert[m]);
  for (i=0; i<min(d,m+1); ++i) C[i]=Csp[i];
}

void sp_hat(int m,int ae[],int v[],int d,int C[])
{
  int m2;			/* No. of vertices */
  int Csp[(max_m+1)];		/* sp_hat(N) weights */
  int invert[max_m+1];		/* Whether the current arc is PQ inverting */
  int type[max_m+1];		/* Surgery type on current arc */
  int isP[2*max_m];		/* Whether current segment is a P */
  int arc_status[max_m];	/* Specific surgery type on current arc */
                                /* 0: regular bridge */
                                /* 1: crossed bridge */
                                /* 2: closed bridge */
  int visited[2*max_m+1];	/* Wilson loop segments visited */
  int seg;			/* segment pointer. Points to the vertex */
				/* before the current segment. */
  int power,sign,anomaly;	/* Power, sign, and anomaly of current cont. */
  int forward;			/* Direction of travel */
  int p;			/* P now? */
  int legit;			/* A legitimate P? */
  int i,j;			/* Temp. var. */

  m2=2*m;
  visited[m2]=1;
  for (i=0; i<=m; ++i) Csp[i]=0;
  for (i=0; i<=m; ++i) invert[i]=0;
  do {
    p=1;
    for (i=0; i<m2; ++i) isP[i]=p^=invert[v[i]];
    legit=1;
    for (i=1; i<m2; ++i) legit=legit && !(invert[v[i]] && isP[i]==isP[ae[i]]);
    if (legit) {
      for (i=0; i<=m; ++i) type[i]=0;
      do {
        for (i=0; i<m2; ++i) visited[i]=0;
        for (i=1; i<m2; ++i) {
          j=v[i];
          if (invert[j]) arc_status[j]=type[j];
          else if (type[j]) arc_status[j]=2;
          else arc_status[j]=(isP[i]!=isP[ae[i]]);
        }
        seg=power=anomaly=0; forward=1; sign=1;
        while (!visited[seg]) {
          ++power;
          while (!visited[seg]) {
            visited[seg]=1;
            switch (arc_status[v[i=(seg+forward) % m2]]) {
            case 0:
              seg=(m2+ae[i]-1+forward) % m2;
              break;
            case 1:
              forward=!forward;
              seg=(m2+ae[i]-1+forward) % m2;
              break;
            case 2:
              seg=(seg+(forward ? 1 : m2-1)) % m2;
              break;
            }
          }
          for (seg=1; seg<m2 && visited[seg]; ++seg) ;
        }

        for (i=0; i<m2; ++i)
            if (!invert[v[i]] && !type[v[i]] && isP[i] && !isP[ae[i]])
               sign= -sign;
        for (i=0; i<m; ++i) if (arc_status[i]==2) {
          sign= -sign;
          ++anomaly;
        }
        --power;
        for (i=0; i<=anomaly; ++i) {
          Csp[power+i]+=sign*binomial[anomaly][i];
          sign*=2;
        }
        for (i=0; type[i]; ++i) type[i]=0;
        type[i]=1;
      } while (!type[m]);
    }
    for (i=0; invert[i]; ++i) invert[i]=0;
    invert[i]=1;
  } while (!invert[m]);
  for (i=0; i<min(d,m+1); ++i) C[i]=Csp[i];
}

	/* WS product computation routine(s) */
	/* Format: product(m,ae,v,n,m_,ws_,d_,C) computes the weight assigned */
	/*         to the diagram (m,ae,v) by the n WSs ws_[0],...,ws_[n] */
	/*         which act on diagrams of m_[0],...,m_[n] loops and */
	/*         produce d_[0],...,d_[n] results each. The result is */
	/*         stored in the vector C whose length is d_[0]*...*d_[n].  */
	/*         The ws_[i]'s are pointers to functions. */

void product(int m,int ae[],int v[],int n,int m_[],
	void (*(ws_[]))(int m,int ae[],int v[],int d,int C[]),
	int d_[],int C[])
{
  int m2;				/* No. of loops */
  int d=1;				/* Total length of C */
  int mask[max_m];			/* `Which arc goes where' mask */
  int m_left[N+1];			/* Arcs left for distribution */
					/* (by types) */
  int ca,cm,*cae,*cv;			/* Current arc, mask, ae, v */
  int nv_[N];				/* No. of vertices treated so */
					/* far for each WS */
  int na_[N];				/* No. of arcs encountered so */
					/* far in each WS */
  int ae_bank[2*max_m],v_bank[2*max_m];
					/* ae and v banks */
  int *(ae_[N]),*(v_[N]);		/* pointer arrays pointing to the */
					/* beginning of the ae's and v's */
  int v_name[2*max_m];			/* New vertex names */
  int *(C_[N]);				/* pointer arrays pointing to the */
					/* C storage area of each WS */
  int *cC;				/* current C storage area */
  int C_bank[max_C_bank];		/* C bank */
  int nC_[N+1];				/* No. of the {C} being treated now */
  int i,j,p;				/* Temporary veriables */

		/* Initializations */
  m2=2*m;
  for (i=0; i<m; ++i) mask[i]=-1;
  for (i=0; i<n; ++i) {
    d*=d_[i];
    m_left[i]=m_[i];
  }
  m_left[n]=1;
  for (i=0; i<d; ++i) C[i]=0;
  v_[0]=cv=v_bank;
  ae_[0]=cae=ae_bank;
  C_[0]=cC=C_bank;
  for (i=1; i<n; ++i) {
    v_[i]=cv+=2*m_[i-1];
    ae_[i]=cae+=2*m_[i-1];
    C_[i]=cC+=d_[i-1];
  }
		/* Producing a mask */
  ca=0;
  while (ca>-1) {
    if (ca<m) {
      cm=mask[ca];
      if (cm>=0) ++(m_left[cm]);
      while (!m_left[++cm] && cm<n);
      if (cm==n) {
        mask[ca]=-1;
        --ca;
      } else {
        --(m_left[cm]);
        mask[ca]=cm;
        ++ca;
      }
    } else {

		/* mask is now produced, now DO SOMETHING!! */
      for (i=0; i<n; ++i) na_[i]=nv_[i]=-1;
      for (i=0; i<m2; ++i) {
        cm=mask[v[i]];
        v_name[i]= ++(nv_[cm]);
        if (ae[i]>i) v_[cm][nv_[cm]]= ++(na_[cm]);
        else {
          v_[cm][nv_[cm]]=v_[cm][v_name[ae[i]]];
          ae_[cm][v_name[ae[i]]]=v_name[i];
          ae_[cm][v_name[i]]=v_name[ae[i]];
        }
      }
      for (i=0; i<n; ++i) {
        (*(ws_[i]))(m_[i],ae_[i],v_[i],d_[i],C_[i]);
      }
      for (i=0; i<=n; ++i) nC_[i]=0;
      do {
        p=1;
        for (i=0; i<n; ++i) p*=C_[i][nC_[i]];
        j=nC_[n-1];
        for (i=n-2; i>=0; --i) j= j*d_[i]+nC_[i];
        C[j]+=p;
        i=0;
        while (i<n && nC_[i]==d_[i]-1) nC_[i++]=0;
        ++(nC_[i]);
      } while (i<n);

		/* End of the mask-producing loop, go to next mask */
      ca=m-2;
      ++(m_left[mask[m-1]]);
      mask[m-1]=-1;
    }
  }
}

		/* Adjoint representations */
		/* Format: adjoint(m,ae,v,d,ws,C) computes the `adjoint' */
		/* to the weight system ws(m,a,v,d) and stores the result */
		/* in C */
void adjoint(int m,int ae[],int v[],int d,
	void (*ws)(int m,int ae[],int v[],int d,int C[]),int C[])
{
  int m2;				/* No. of vertices */
  int left[2*max_m+1];			/* left turn signal */
  int ae_ws[2*max_m+2],v_ws[2*max_m+2];
					/* ae and v in ws */
  int beginning[max_m];		/* Beginning of arc #i */
  int C_ws[max_m+2];			/* ws's C */
  int cv,cl,cr;				/* Current vertex, left, right */
  int lQ;				/* left now? */
  int sign;				/* Sign of current term */
  int i,j;				/* Temporary variable */

  m2=2*m;
  for (i=0; i<d; ++i) C[i]=0;
  for (i=0; i<=m2; ++i) left[i]=0;
  do {
    sign=1;
    cl=0; cr=m2;
    for (i=0; i<m2; ++i) {
      cv=((lQ=left[i]) ? cl : cr);
      if (i < ae[i]) beginning[v_ws[cv]=v[i]]=cv;
      else {
        j=ae_ws[cv]=beginning[v_ws[cv]=v[i]];
        ae_ws[j]=cv;
      }
      if (lQ) ++cl;
      else {
        --cr;
        sign=-sign;
      }
    }
    v_ws[cl]=v_ws[m2+1]=m;
    ae_ws[cl]=m2+1; ae_ws[m2+1]=cl;
    (*ws)(m+1,ae_ws,v_ws,d,C_ws);
    for (i=0; i<d; ++i) C[i]+=sign*C_ws[i];
    for (i=0; left[i]; ++i) left[i]=0;
    left[i]=1;
  } while (!left[m2]);
}
// ws.h"
#define Lie_print               0       /* Print Lie solutions? */
#define	print_all_equations	1	/* Print all equations? */
					/* 0 - no */
					/* 1 - yes, all */
					/* 2 - yes, basis only */

void print_equation(struct EQUATION equation);
void print_polynom(char var, int deg, int *coeffs);

void print_diagram(int m,int count)
{
  int v[2*M];			/* v */
  int C[max_m+1];		/* C storage */
  int m2;			/* No. of vertices */
  int m_[2]; 
  void fuck_you(), ((*(ws_[2]))(int , int [], int [], int , int []));
				/* fuck_you is a proper compiler directive */
				/* here. */
  int d_[2];
				/* `product' parameters */
  int i;			/* temporary variable */
  struct EQUATION eqn;          /* Temporary storage for */
                                /* compressed acc */
  int ae[2*M];			/* ae */
  m2=2*m;

  uncompress_ae(m,diags[count],ae);
  make_v(m,ae,v);
  printf("\n\n");
  print_d(m,count); printf("="); printCD(v);

  if (Lie_print) {
    printf("\nW[su,N_,def][");
    print_d(m,count);
    printf("]:=");
    su(m,ae,v,(m+3)/2,C);
    print_polynom('N',m/2,C);

    printf("\nW[so,N_,def][");
    print_d(m,count);
    printf("]:=");
    so(m,ae,v,m+1,C);
    print_polynom('N',m,C);

    printf("\nW[sp,N_,def][");
    print_d(m,count);
    printf("]:=");
    sp(m,ae,v,m+1,C);
    print_polynom('N',m,C);

    printf("\nWH[su,N_,def][");
    print_d(m,count);
    printf("]:=");
    su_hat(m,ae,v,(m+3)/2,C);
    print_polynom('N',m/2,C);

    printf("\nWH[so,N_,def][");
    print_d(m,count);
    printf("]:=");
    so_hat(m,ae,v,m+1,C);
    print_polynom('N',m,C);

    printf("\nWH[sp,N_,def][");
    print_d(m,count);
    printf("]:=");
    sp_hat(m,ae,v,m+1,C);
    print_polynom('N',m,C);
  }

  if (solveQ) {
    for (k=0; k<ndiags; ++k) acc.coeff[k]=0;
    k=ndiags;
    for (i=0; i<afs_depth; ++i) {
      k>>=afs_step;
      for (j=0; j<=k; ++j) afs[i][j]=0;
    }
    acc.len=1; afs_update(count,1);
    acc.coeff[count]=1;
    to_basis(ndiags,ndiags); eqn=compress(acc);
    if (framed) printf("\n\nToBasisA[");
    else printf("\n\nToBasisAr[");
    print_d(m,count);
    printf("]:=");
    print_equation(eqn);
  }
}

                /* Printing the solution */
void print_solution()
{
  int i,j,k,deg;			/* Temporary variables */
  int m2;               	        /* No. of vertices */
  int ae[2*M];                          /* The other end of the arc */

  m2=2*M;
  printf("\n");
  if (print_all_equations)
     for (count=0; count<ndiags; ++count) 
	 if (print_all_equations==1 || !matrix[count].len)
	    print_diagram(M,count);

                /* Computing the number of degrees of freedom */
  if (solveQ) {
    deg=0;
    if (framed) printf("\nBasisA[%d]={",M);
    else printf("\nBasisAr[%d]={",M);
    for (i=0; i<ndiags; ++i) if (!matrix[i].len) {
      if (deg!=0) printf(",");
      if (deg%5==4) printf("\\\n	");
      uncompress_ae(M,diags[i],ae);
      make_v(M,ae,v);
      print_d(M,i);
      ++deg;
    }
    printf("}");
    printf("\n(*	Total nunmber of diagrams is: %6d.	*)\n",ndiags);
    printf("(*	Number of degrees of freedom is:%4d.	*)\n",deg);
  }
}

void print_equation(struct EQUATION equation)
{
  int i;

  if (equation.len==0) printf("0");
  else {
    for (i=0; i<min(ndiags,equation.len); ++i) {
      if (i%5 == 3) printf("\\\n	");
      if (i>0 && (!equation.coeff[i]) > 0) printf("+");
      if ((abs(!equation.coeff[i]))!=1) printf("%d",!equation.coeff[i]);
      else if ((!equation.coeff[i])==-1) printf("-");
      print_d(M,equation.at[i]);
    }
  }
}

void print_monom(char var, int power, int coeff)
{
  if (power==0) printf("%d",coeff);
  else {
    if (coeff==1) ;
    else if (coeff==-1) printf("-");
    else printf("%d",coeff);
    printf("%c",var);
    if (power>1) printf("^%d",power);
  }
}

void print_polynom(char var, int deg, int *coeffs)
{
  int i;

  i=deg;
  while (i>=0 && coeffs[i]==0) --i;
  if (i>=0) print_monom(var,i,coeffs[i]); else printf("0");
  --i;
  while (i>=0) {
    while (i>=0 && coeffs[i]==0) --i;
    if (i>=0) {
      if (coeffs[i]>0) printf("+");
      print_monom(var,i,coeffs[i]);
      --i;
    }
  }
}

void printCD(int v[2*M])
{
  int i;

  printf("CD[%d",v[0]+1);
  for (i=1; i<2*M; ++i) printf(",%d",v[i]+1);
  printf("]");
}

void print_d(int m,int count)
{
  printf("d");
  if (!framed) printf("r");
  printf("[%d,%d]",m,count+1);
}

int statisticsQ()
{
  int c,diff;		/* Clock reading */
  static int last_c=0;	/* internal rep */

  c=clock();
  diff=c-last_c;
  if (diff>1000000*statistics_interval) {
    ++rep;
    printf("(*t=%d ",rep*statistics_interval);
    last_c=c;
    return(rep);
  } else return(0);
}
