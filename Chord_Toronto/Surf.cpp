#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>

typedef int COEFF;	// The ground ring.

		// The Chinese Character container class
class CC {
private:
  int m;		// The overall degree
  int e;		// Number of univalent vertices
  int *oea;		// The other end of end #$
  int *arca;		// The arc connected to end #$
  int *uva;		// The end of univalent vertex #$
  int *aba;		// The beginning of iarc #$
  int *aea;		// The end of iarc #$
  void update();	// Update after a change in oe

public:
  CC(int m, int e);			// Constructor
  CC(const CC &cc);			// The copy constructor
  ~CC();				// Destructor
  int mQ() const {return m;}
  int eQ() const {return e;}
  int oe(int i) const {return oea[i];}	// The other end of end #$
  int arc(int i) const {return arca[i];}// The arc connected to end #$
  int uv(int i) const {return uva[i];}	// The end of univalent vertex #$
  int ab(int i) const {return aba[i];}	// The beginning of iarc #$
  int ae(int i) const {return aea[i];}	// The end of iarc #$
  void toH(int at_arc);			// change I->H at iarc at_arc
  void flip(int at_vertex);		// flip the orientation of a vertex
};

ostream& operator << (ostream& s, CC cc);	// printing out a CC
int No2UVs(const CC& cc);	// Checks that no 2 UVs are one the save vertex
int No1loops(const CC& cc);	// No 1-arc loops
int No2loops(const CC& cc);	// No 2-arc loops
int No3loops(const CC& cc);	// No 3-arc loops

		// The marked surface container classes
class CSURF;
class SURF {
private:
  int o;		// Orientability of our SURF
  int max_mark;		// Current maximal marking number; default -1.
  int *b;		// b[0] unmarked boundaries,
			// b[1] boundaries marked once,...
			// b[max_mark] boundaries marked "max_mark" times.

public:
  SURF(int oriented);			// orientation-dependant constructor
  SURF(const SURF &surf);		// The copy constructor
  SURF(CSURF csurf);			// The SURF part of a CSURF.
  ~SURF();				// Destructor
  int orientedQ() {return o;}		// orientability inquiry
  int bot(int k) {return b[k];}		// # of k-marked boundaries
  int max_marking() {return max_mark;}	// returns max_mark
  int add_bot(int k);			// add one bot k
};

ostream& operator << (ostream& s, SURF surf);	// printing out a SURF
int compare(SURF surf1, SURF surf2);	// Lex. ordering of surf1, surf2:
					// surf1 < surf2 : 1; =:0; >:-1.

class CSURF : public SURF {
public:
  COEFF coeff;				// The coefficient
  CSURF(int oriented) :SURF(oriented), coeff(1) {}
			// orientation-dependant constructor
  //~CSURF() {SURF::~SURF();}		// Destructor
  CSURF& operator *= (COEFF factor) {coeff*=factor; return *this;}
					// Multiplying a CSURF by a const.
};

ostream& operator << (ostream& s, CSURF csurf);	// printing out a CSURF

class SURF_COMB {
private:
  int sz;			// The total length of the surface combination
  CSURF *middle;		// The "middle" CSURF.
  SURF_COMB *left,*right;	// Right and left branches

public:
  SURF_COMB();				// Constructs the zero combination
  ~SURF_COMB();				// Destructor
  int size() {return sz;}		// returns "sz"
  SURF_COMB& operator += (CSURF csurf);	// Adds in a CSURF
  CSURF *(take_left());		// Take (and remove) the left most branch
  CSURF *(take_right());	// Take (and remove) the right most branch
  friend ostream& operator << (ostream& s, SURF_COMB sc);
				// printing out a SURF_COMB
};

                // Cyclically next & previous ends
inline int next_end(int i) {if (i%3<2) return (i+1); else return (i-2);}  
inline int prev_end(int i) {if (i%3>0) return (i-1); else return (i+2);}

		// The Phi functions
CSURF tau(const CC &cc, int crossed[]);	// Thickening a CC into a SURF
SURF_COMB& Phi(CC cc);			// The Phi map.
// surf.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

main() {
  int at_arc;		// Do toH at at_arc.
  int at_vertex;	// flip vertex # at_vertex
  int m,e;

  cout << "(* Enter m: "; cin >> m;	// The overall degree
  cout << "m=" << m << " *)\n";
  cout << "(* Enter e: "; cin >> e;	// Number of univalent vertices
  cout << "e=" << e << " *)\n";
  cout << "(* Enter number of CC's to construct: ";
  int n;
  cin >> n;
  cout << "n=" << n << " *)\n";
  cout << "(* Enter number of random moves after each Phi computation: ";
  int k;
  cin >> k;
  cout << "k=" << k << " *)\n";

  cout << "\nl[" << m << "," << e << "]={};\n\n";
  CC cc(m,e);		// The Chinese Character presently being considered
  for (int i=0; i<n; ++i) {
    if (No2UVs(cc) && No1loops(cc) &&
                      (i%20==0 || (No2loops(cc) && No3loops(cc)))) {
      cout << "(* " << i << " " << cc << " *)\n";
      cout.flush();
      SURF_COMB res=Phi(cc);
      if (res.size()==0) cout << 0 << ";\n\n";
      else cout << "AppendTo[l[" << m << "," << e << "],\n" << res << "];\n\n";
      cout.flush();
    }
    for (int j=0; j<k; ++j) {
      at_vertex=((int) random())%(2*m-e);
      cc.flip(at_vertex);
      at_arc=(int) random();
      at_arc%=(3*m-2*e);
      while (cc.ab(at_arc)/3 == cc.ae(at_arc)/3) at_arc=(at_arc+1)%(3*m-2*e);
      //cout << "toH(" << at_arc << "): ";
      cc.toH(at_arc);
    }
  }
}


// surf.h"
#include <iostream>
#include <fstream>
#include <iomanip>

CC::CC(int m_in, int e_in) 	// The CC constructor
{
  m=m_in; e=e_in;
  oea=new int[6*m-3*e]; arca=new int[6*m-3*e];
  uva=new int[e];
  aba=new int[3*m-2*e]; aea=new int[3*m-2*e];
  int i;
  for (i=1; i<2*m-e; ++i) {oea[3*i-1]=3*i; oea[3*i]=3*i-1;}
  oea[6*m-3*e-1]=0; oea[0]=6*m-3*e-1;
  for (i=0; i<e; ++i) oea[3*i+1]=3*i+1;
  for (i=e; i<2*m-e; i+=2) {oea[3*i+1]=3*i+4; oea[3*i+4]=3*i+1;}
  update();
}

CC::CC(const CC &cc)	// The copy constructor
{
  m=cc.m; e=cc.e;
  oea=new int[6*m-3*e]; arca=new int[6*m-3*e];
  uva=new int[e];
  aba=new int[3*m-2*e]; aea=new int[3*m-2*e];
  for (int i=0; i<6*m-3*e; ++i) oea[i]=cc.oea[i];
  update();
}

CC::~CC()		// Destructor
{
  delete oea; delete arca; delete uva; delete aba; delete aea;
}

void CC::update() // Update after a change in oe
{
  int i;
  for (i=0; i<6*m-3*e; ++i) arca[i]=-1;
  int iarc_enc=0;	// # of iarcs encountered so far
  int uv_enc=0;		// # of univalent vertices encountered so far
  for (i=0; i<6*m-3*e; ++i) if (arca[i]==-1) {
    if (i==oea[i]) {
      uva[uv_enc]=i;
      arca[i]=3*m-2*e+uv_enc;
      ++uv_enc;
    } else {
      arca[i]=arca[oea[i]]=iarc_enc;
      aba[iarc_enc]=i;
      aea[iarc_enc]=oea[i];
      ++iarc_enc;
    }
  }
}

/* ostream& operator << (ostream& s, CC cc)	// printing out a CC
{
  s << "CC[";
  int i;
  for (i=0; i<2*m-e;) {
    s << "OV[" << cc.arc(3*i) << "," << cc.arc(3*i+1) << ","
         << cc.arc(3*i+2) << "]";
    if ((++i)<2*m-e) s << ",";
  }
  for (i=0; i<e; ++i) s << ",UV[" << cc.arc(cc.uv(i)) << "]";
  return s << "]";
} */

ostream& operator << (ostream& s, CC cc)	// printing out a CC
{
  s << "{";
  int i;
  int m=cc.mQ(); int e=cc.eQ();
  for (i=0; i<2*m-e;) {
    s << cc.arc(3*i) << " " << cc.arc(3*i+1) << " " << cc.arc(3*i+2);
    if ((++i)<2*m-e) s << "; ";
  }
  if (e>0) s << " ; ";
  for (i=0; i<e; ++i) s << " " << cc.arc(cc.uv(i));
  return s << "}";
}

void CC::toH(int at_arc)		// change I->H at arc at_arc
{
  int ne=next_end(aba[at_arc]);		// NE corner of I/H
  int sw=prev_end(aea[at_arc]);		// SW corner of I/H
  if (oea[sw]==ne) return;
  oea[oea[ne]]=sw;
  oea[oea[sw]]=ne;
  int t=oea[ne]; oea[ne]=oea[sw]; oea[sw]=t;
  update();
}

void CC::flip(int at)			// flip the orientation of a vertex
{
  at*=3;
  if (oea[at+1]==at+2 || oea[at+1]==at+1 && oea[at+2]==at+2) return;
  if (oea[at+1]==at+1) {
    oea[at+1]=oea[at+2];
    oea[oea[at+1]]=at+1;
    oea[at+2]=at+2;
    update(); return;
  }
  if (oea[at+2]==at+2) {
    oea[at+2]=oea[at+1];
    oea[oea[at+2]]=at+2;
    oea[at+1]=at+1;
    update(); return;
  }
  int t=oea[at+1]; oea[at+1]=oea[at+2]; oea[at+2]=t;
  oea[oea[at+1]]=at+1;
  oea[oea[at+2]]=at+2;
  update();
}

int No2UVs(const CC& cc)	// Checks that no 2 UVs are one the save vertex
{
  int m=cc.mQ(); int e=cc.eQ();
  for (int i=0; i<2*m-e; ++i) {
    int uvs=0;
    for (int j=0; j<3; ++j) if (3*i+j==cc.oe(3*i+j)) ++uvs;
    if (uvs>=2) return 0;
  }
  return 1;
}

int No1loops(const CC& cc)	// No 1-arc loops
{
  int m=cc.mQ(); int e=cc.eQ();
  for (int i=0; i<3*m-2*e; ++i) if (cc.ab(i)/3==cc.ae(i)/3) return 0;
  return 1;
}

int No2loops(const CC& cc)	// No 2-arc loops    
{
  int m=cc.mQ(); int e=cc.eQ(); int at;
  for (int i=0; i<6*m-3*e; ++i) if (i/3!=(at=cc.oe(i)))
      if (i/3==cc.oe(next_end(at))/3 || i/3==cc.oe(prev_end(at))/3)
         return 0;
  return 1;
}

int No3loops(const CC& cc)	// No 3-arc loops
{
  int m=cc.mQ(); int e=cc.eQ(); int atp,atq,atr;
  for (int i=0; i<6*m-3*e; ++i) if (i!=(atp=cc.oe(i))) for (int p=0; p<2; ++p)
  {
    atp=(p ? next_end(atp) : prev_end(atp));
    if (atp!=(atq=cc.oe(atp))) for (int q=0; q<2; ++q) {
      atq=(q ? next_end(atq) : prev_end(atq));
      if (atq!=(atr=cc.oe(atq))) if (atr/3==i/3) return 0;
    }
  }
  return 1;
}
// surf.h"

SURF::SURF(int oriented)	// The SURF orientation-dependant constructor
{
  max_mark=-1;
  o=oriented;
}

SURF::SURF(const SURF &surf)	// The copy constructor
{
  o=surf.o;
  if (surf.max_mark>=0) {
    b=new int[1+(max_mark=surf.max_mark)];
    for (int i=0; i<=max_mark; ++i) b[i]=surf.b[i];
  } else max_mark=-1;
}

SURF::SURF(CSURF csurf)		// The SURF part of a CSURF.
				// ugly implementation.
{
  o=csurf.orientedQ();
  max_mark=csurf.max_marking();
  if (max_mark>=0) {
    b=new int[max_mark+1];
    for (int i=0; i<=max_mark; ++i) b[i]=csurf.bot(i);
  }
}

SURF::~SURF()		// Destructor
{
  if (max_mark>=0) delete b;
}

int SURF::add_bot(int k)	// add one bot k
{
  if (k<=max_mark) return ++(b[k]);
  if (max_mark>=0) {
    int *bt=b;
    b=new int[k+1];
    for (int i=0; i<=max_mark; ++i) b[i]=bt[i];
    delete bt;
  } else b=new int[k+1];
  for (int i=max_mark+1; i<k; ++i) b[i]=0;
  return b[max_mark=k]=1;
}

//ostream& operator << (ostream& s, SURF surf)	// printing out a SURF
//{
  //s << "MarkedSurface[";
  //if (surf.orientedQ()) s << "Orientable"; else s << "NonOrientable";
  //for (int i=0; i<=surf.max_marking(); ++i) s << "," << surf.bot(i);
  //return s << "]";
//}

ostream& operator << (ostream& s, SURF surf)	// printing out a SURF
{
  if (surf.orientedQ()) s << "OY"; else s << "ON";
  s << "[";
  s << surf.bot(0);
  for (int i=1; i<=surf.max_marking(); ++i) s << "," << surf.bot(i);
  return s << "]";
}

int compare(SURF surf1, SURF surf2)	// Lex. ordering of surf1, surf2:
					// surf1 < surf2 : 1; =:0; >:-1.
{
  if (surf1.orientedQ()!=surf2.orientedQ()) if (surf1.orientedQ()) return -1;
     else return 1;
  if (surf1.max_marking()!=surf2.max_marking())
     if (surf1.max_marking()<surf2.max_marking()) return 1; else return -1;
  for (int i=0; i<=surf1.max_marking(); ++i) if (surf1.bot(i)!=surf2.bot(i))
                           if (surf1.bot(i)<surf2.bot(i)) return 1;
                           else return -1;
  return 0;
}

ostream& operator << (ostream& s, CSURF csurf)	// printing out a CSURF
{
  return s << "(" << csurf.coeff << ")" << ((SURF) csurf);
}

SURF_COMB::SURF_COMB()		// Constructs the zero combination
{
  sz=0;
  left=right=0;
  middle=0;
}

SURF_COMB::~SURF_COMB()		// Destructor
{
  if (left) {delete left; left=0;}
  if (middle) {delete middle; middle=0;}
  if (right) {delete right; right=0;}
}

SURF_COMB& SURF_COMB::operator += (CSURF csurf)	// Adds in a CSURF
{
  if (csurf.coeff==0) return *this;
  if (sz==0) {
    sz=1;
    middle=new CSURF(csurf);
    left=new SURF_COMB;
    right=new SURF_COMB;
    return *this;
  }
  int c=compare(*middle,csurf);
  if (c!=0) {
    if (c==1) (*right)+=csurf; else (*left)+=csurf;
    sz=1+left->size()+right->size();
  } else {
    middle->coeff+=csurf.coeff;
    if (middle->coeff==0) {
      sz--;
      delete middle; middle=0;
      if (left->size()>right->size()) {middle=left->take_right(); 
			return *this;}
      if (left->size()<=right->size() && right->size()>0)
         {middle=right->take_left();
			 return *this;}
      delete left; delete right;
      left=right=0;
    }
  }
}

		// Take (and remove) the left most branch
CSURF *(SURF_COMB::take_left()) {
  --sz;
  if (sz==0) {delete left; delete right; CSURF *t=middle;
    left=right=0; middle=0; return t;}
  if (left->size()>0) return left->take_left();
  delete left;
  left=right->left; CSURF *t=middle; middle=right->middle; right=right->right;
  return t;
}

		// Take (and remove) the right most branch
CSURF *(SURF_COMB::take_right()) {
  --sz;
  if (sz==0) {delete left; delete right; CSURF *t=middle;
    left=right=0; middle=0; return t;}
  if (right->size()>0) return right->take_right();
  delete right;
  right=left->right; CSURF *t=middle; middle=left->middle; left=left->left;
  return t;
}

ostream& operator << (ostream& s, SURF_COMB sc)	// printing out a SURF_COMB
{
  if (sc.size()==0) return s;
  if (sc.left!=0) s << *sc.left;
  s << "+" << *sc.middle;
  if (sc.right!=0) s << *sc.right;
  return s;
}
// surf.h"

	// Recursively verify orientability, beginning at vertex #from.
int verify_orientability(const CC &cc, int *crossed, int from, int *up_side);

CSURF tau(const CC &cc, int crossed[])	// Thickening a CC into a CSURF
					// Warning!! bug if cc=Y
{
  int m=cc.mQ(); int e=cc.eQ();
  int *up_side=new int[2*m-e];	// Which side is up in vertex #$.
  int *visited=new int[6*m-3*e];
			// Whether armpit #$ has been traveled through
  int i;

  up_side[0]=1;
  for (i=1; i<2*m-e; ++i) up_side[i]=-1;
  CSURF csurf(verify_orientability(cc,crossed,0,up_side));
			// Creating the surface.
  for (i=0; i<6*m-3*e; ++i) visited[i]=0;

  for (int start=0; start<6*m-3*e; ++start) if (!visited[start]) {
    int markings=0;	// markings encountered so far
    int direction=(csurf.orientedQ() ? up_side[start/3] : 1);
			// current direction of travel
    int at=start;
    do {
      visited[at]=1;
      if (direction) at=next_end(at);
      if (at==cc.oe(at)) {
        ++markings;
        if (!direction) csurf.coeff*=-1;
        if (!direction) at=prev_end(at);
      } else {
        direction^=crossed[cc.arc(at)];
        at=cc.oe(at);
        if (!direction) at=prev_end(at);
      }
    } while (!visited[at]);
    csurf.add_bot(markings);
    if (!csurf.orientedQ() && markings%2) csurf.coeff=0;
  }
  if (e%2) csurf.coeff=0;
  delete up_side; delete visited;
  return csurf;
}

	// Recursively verify orientability, beginning at vertex #from.
int verify_orientability(const CC &cc, int *crossed, int from, int *up_side)
{
  for (int i=3*from; i<3*from+3; ++i) if (up_side[cc.oe(i)/3]==-1) {
    up_side[cc.oe(i)/3]=crossed[cc.arc(i)]^up_side[from];
    if (!verify_orientability(cc,crossed,cc.oe(i)/3,up_side)) return 0;
  } else if (i!=cc.oe(i) &&
             crossed[cc.arc(i)]^up_side[from]^up_side[cc.oe(i)/3]) return 0;
  return 1;
}

SURF_COMB& Phi(CC cc)		// The Phi map.
{
  int m=cc.mQ(); int e=cc.eQ();
  int *crossed=new int[3*m-2*e+1];	// Whether $'th arc is crossed or not.

  SURF_COMB *res=new SURF_COMB;
  int i;
  for (i=0; i<=3*m-2*e; ++i) crossed[i]=0;
  COEFF sign=1;
  do {
    (*res)+=(tau(cc,crossed)*=sign);
    for (i=0; crossed[i]; ++i) {sign*=-1; crossed[i]=0;}
    sign*=-1; crossed[i]=1;
  } while (i<3*m-2*e);
  delete crossed;
  return *res;
}
