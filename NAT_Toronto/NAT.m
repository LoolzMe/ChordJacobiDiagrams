(* (c) 12/20/93, Dror Bar-Natan. *)
BeginPackage["NAT`",{"LinAlg`"}]

$RecursionLimit = Max[65536,$RecursionLimit]

t::usage = "t[i,j] is a horizontal chord connecting the i'th strand with \
	the j'th strand."
p::usage = "p[c1,c2,...] represents the product of the horizontal \
	chords c1,c2,..."
b::usage = "b[x,y] represents the Lie bracket of x and y."

ModDegree::usage = "ModDegree[m,expr] computes expr ignoring everything \
	whose degree is m or higher"

Unprotect[NonCommutativeMultiply]
ClearAll[NonCommutativeMultiply]
NonCommutativeMultiply::usage = "d1**d2 computes the product of d1 and d2 in AP"

APPower::usage = "APPower[x,n] is x raised to the nth power using **"

d::usage = "d[n][expr] computes the nth differential of expr. d[n,i][expr] \
	is the (n,i) coface of expr. d and OX interact correctly; \
	expressions like d~OX~1 are defined correctly."

OX::usage = "The formal tensor product symbol."

APReduce::usage = "APReduce[expr] reduces all p-products in expr modulo the \
	Four Term relation to a cannonical form"

Act::usage = "Act[spec][expr] acts on all the t[i,j]'s in expr according \
	to the specification spec"

AP2CD::usage = "AP2CD[ord,expr] takes an element of the algebra AP_n and \
	closes it to a linear combination of chord diagrams in A^r by \
	traveling around the strands of AP_n in the order specified by \
	ord, going up at first and then down, up, down, ..."

CD::usage ="CD[chords] represents a chord diagram"

Invert::usage = "Invert[1+expr] inverts using a formal power series"

FormalLog::usage = "FormalLog[1+expr] computes the logarithm of (1+expr) \
	using a power series expansion with **."

T::usage = "T[n] is the transpose operation of AP_n"

NormSquared::usage = "NormSquared[expr] is the norm-squared of an element \
	of AP_3"

CDReduce::usage = "CDReduce[expr] will rewrite all chord diagrams in expr \
	in terms of a basis of A^r"

AllCDs::usage = "AllCDs[m] is a List of all the chord diagrams in degree m"

HarrisonBasis::usage = "HarrisonBasis[n,m] returns a basis to the degree m \
	piece of C_Harr^n."

TotalBracket::usage = "TotalBracket[n] is the n'th total super-bracketing \
	operator, whose image is HarrisonBasis[n,m]"

CBasis::usage = "CBasis[n,m] returns a basis to the degree m \
	piece of C^n"

LieExpand::usage = "LieExpand[expr] replaces all Lie brackects in expr \
	by commutators."

LieBasis::usage = "LieBasis[m,gens] returns a basis of the space of Lie \
	polnomials of degree m in the generators gens."

FullLieBasis::usage = "FullLieBasis[m,gens] returns a basis of the space \
	of Lie polnomials of degrees up to m in the generators gens."

EvenLieBasis::usage = "EvenLieBasis[m,gens] returns a basis of the space \
	of Lie polnomials of even degrees up to m in the generators gens."

WHOMFLY::usage = "WHOMFLY[N,h][cd] computes the weight of the \
	chord diagram cd in the weight system corresponding to the \
	HOMFLY polynomial. Here h is log(q) and N is the N of SL(N)"

Begin["`private`"]

LieExpand[expr_]:=Expand[expr //. b[x_,y_] :> x**y-y**x]

LieBasis[1,gens_] := gens
LieBasis[m_,gens_] := LieBasis[m,gens] = Module[{lst},
  lst=Flatten[Outer[b,gens,LieBasis[m-1,gens]]];
  lst[[SpanningIndices[LieExpand /@ lst]]]
]

FullLieBasis[m_,gens_] := Join @@ Table[LieBasis[i,gens],{i,m}]

EvenLieBasis[m_,gens_] := Join @@ Table[LieBasis[i,gens],{i,2,m,2}]

If[Head[mstack]=!=List,mstack={Infinity}]
currentm=Last[mstack]

SetAttributes[ModDegree,HoldRest]
ModDegree[m_,expr_] := Module[{res},
  AppendTo[mstack,currentm=m];
  res=expr;
  mstack=Drop[mstack,-1];
  currentm=Last[mstack];
  res
]

If[Head[DimCDr[1]]===DimCDr,Get["CDReduceData.m"]]

CDReduce[expr_] := Expand[
  expr /. cd_CD :> Module[
    {lcd,f,i},
    If[Length[cd]>14,cd,
      If[Or @@ MapThread[Equal,{List @@ cd,RotateLeft[List @@ cd]}],0,
        lcd=cd;
        While[Head[ToBasisAr[lcd]]===ToBasisAr,
          f=First[lcd=RotateRight[lcd]];
          lcd=(lcd /. f->0) /. i_ :> i+1 /; i<f
        ];
        ToBasisAr[lcd]
      ]
    ]
  ]
]

AllCDs[m_] := If[m>7,
  Print["Too high degree"]; {},
  Table[dr[m,i], {i,1,DimCDr[m]}]
]

x_CD ** y_CD := If[Length[x]+Length[y] < 2 currentm,
  CDReduce[x ~Join~ (Length[x]/2+# & /@ y)],
  0
]

APPower[x_,n_] := NonCommutativeMultiply @@ Table[x,{n}]

NonCommutativeMultiply[]=1
NonCommutativeMultiply[x_]:=x

c_?NumberQ**x_ := c*x	(* multiplication by a scalar *)
x_**c_?NumberQ := c*x
		(* bilinearity of NonCommutativeMultiply: *)
(c_?NumberQ*x_)**y_ := Expand[c*(x**y)]
x_**(c_?NumberQ*y_) := Expand[c*(x**y)]
(x_+y_)**z_ := x**z + y**z
x_**(y_+z_) := x**y + x**z
		(* multiplying two products: *)
p1_p**p2_p := If[Length[p1]+Length[p2] < currentm,APReduce[p1 ~Join~ p2],0]
t1_t**p2_p := If[Length[p2]+1 < currentm,APReduce[p[t1] ~Join~ p2],0]
p1_p**t2_t := If[Length[p1]+1 < currentm,APReduce[p1 ~Join~ p[t2]],0]
t1_t**t2_t := If[2 < currentm,APReduce[p[t1,t2]],0]
		(* multiple multiplication: *)
NonCommutativeMultiply[x__] := Fold[
  NonCommutativeMultiply,
  First[{x}],
  Drop[{x},1]
] /; Length[{x}] > 2
(*	NonCommutativeMultiply[lft___,mid_Plus,rgt___] := Apply[
	  NonCommutativeMultiply,
	  Distribute[seq[lft,mid,rgt]],
	  {1}
	]	*)

APRules={
  p[lft___,t[i_,j_],t[k_,l_],rgt___] :>
    Which[
      i==k, p[lft,t[k,l],t[i,j],rgt] + p[lft,t[l,j],t[i,j],rgt] -
            p[lft,t[i,j],t[l,j],rgt],
      i==l, p[lft,t[k,l],t[i,j],rgt] + p[lft,t[k,j],t[i,j],rgt] -
            p[lft,t[i,j],t[k,j],rgt],
      True, p[lft,t[k,l],t[i,j],rgt]
    ] /; (j>l)
}
APReduce[expr_]:=FixedPoint[
  (Expand[# /. APRules])&,
  Expand[expr /. t[i_,j_] :> t[j,i] /; i>j]
]

Act[spec___][expr_] := Module[{mspec},
  mspec=Replace[#,c_Integer:>{c}]& /@ {spec};
  APReduce[
    expr /. t[i_,j_] :> Plus @@ Plus @@ Outer[t,mspec[[i]],mspec[[j]]]
         /. p[lft___,0,rgt___] -> 0
         /. x_p :> Distribute[x]
  ]
]

AP2CD[ord_List,c_?NumberQ]=c
AP2CD[ord_List,x_ + y_] := AP2CD[ord,x] + AP2CD[ord,y]
AP2CD[ord_List,c_?NumberQ x_] := Expand[c AP2CD[ord,x]]
AP2CD[ord_List,x_t] := AP2CD[ord,p[x]]

AP2CD[ord_List,x_p] := CDReduce[Module[
  {tc=Table[0,{Length[ord]}],X=List @@ x,cd,i,l,tt,k,j,sign=1},
  X=X /. Thread[Rule[ord,Range[Length[ord]]]];
  X=X /. i_?NumberQ :> {i,++tc[[i]]};
  tt=FoldList[Plus,0,tc];
  X=Sort[X /. {
    {i_?OddQ, l_} :> tt[[i]]+l,
    {i_?EvenQ,l_} :> (sign*=-1; tt[[i]]+tc[[i]]-l+1)
  }];
  k=0;
  X=X /. t[i_,j_] :> t[{i,++k},{j,k}];
  X=Last /@ Sort[Flatten[t @@ X]];
  sign*(CD @@ X)
]]

Invert[1]=1
Invert[1+expr_] := Module[{term,s},
  For[
    term=Expand[-expr] ; s=1,
    !term===0,
    term=term**(-expr),
    s+=term
  ];
  s
]

FormalLog[1]=0
FormalLog[1+expr_] := Module[{term,s,i},
  For[
    i=1; term=Expand[expr] ; s=0,
    term=!=0,
    ++i; term=term**(-expr),
    s+=Expand[term/i]
  ];
  s
]

T[n_][expr_] := OX[T[n]][expr]

d[n_][expr_] := OX[d[n]][expr]

NormSquared[expr_] := expr**T[3][expr]

SetAttributes[OX,Flat]
(OX[l___] /; MemberQ[OX[l],d[_]])[expr_] := Through[
  Distribute[
    OX[l] /. d[n_] :> Sum[(-1)^k * d[n,k],{k,0,n+1}]
  ][expr]
]
OX[lft___,-mid_,rgt___][expr_] := -OX[lft,mid,rgt][expr]
OX[l___][expr_] := Module[
  {
    shifts=FoldList[Plus,0,Drop[{l},-1] /. {T[n_] :> n, d[n_,_] :> n+1}],
    actions={l} /. {
      T[n_] :> Range[n,1,-1],
      d[n_,0] :> Range[2,n+1],
      d[n_,i_] :> Range[n] /; i==n+1,
      d[n_,i_] :> Range[i-1]~Join~{{i,i+1}}~Join~Range[i+2,n+1] /; 0<i<=n,
      n_Integer :> Range[n]
    }
  },
  (Act @@ (Join @@ (actions+shifts)))[expr]
]

HarrisonBasis[n_,m_] := HarrisonBasis[n,m] = LinearSpan[
  TotalBracket[n] /@ (LieExpand[LieBasis[m,CBasis[n,1]]]~Join~CBasis[n,m])
]

TotalBracket[n_][expr_]:=Module[{i,tb=expr},
  Do[
    tb+=Expand[(-1)^i (Act @@ (RotateLeft[Range[i]]~Join~Range[i+1,n]))[tb]],
    {i,2,n}
  ];
  tb
]

CBasis[n_,0] = {1}
CBasis[n_,1] := CBasis[n,1] = Flatten[Table[t[i,j],{i,n-1},{j,i+1,n}]]
CBasis[n_,m_] := CBasis[n,m] = ModDegree[m+1,Module[
  {bas,k,i,j,l},
  bas=CBasis[n,m-1];
  Sort[Flatten[Table[
    l=Last[bas[[k]]];
    If[!IntegerQ[l],l=Last[l]];
    Table[bas[[k]]**t[i,j],{j,l,n},{i,j-1}],
    {k,Length[bas]}
  ]]]
]]

CD2CDP[cd_CD] := Module[
  {m=Length[cd]/2,i,beg,end},
  beg=Table[Position[cd,i],{i,1,m}];
  end=Reverse /@ beg;
  CDP @@ (Range[2m] /. Thread[Rule[Flatten[beg],Flatten[end]]])
]

WHOMFLY[N_,h_][c_?NumberQ] := c*N
WHOMFLY[N_,h_][c_?NumberQ*cd_CD] := c*WHOMFLY[N,h][cd]
WHOMFLY[N_,h_][x_+y_]:=WHOMFLY[N,h][x]+WHOMFLY[N,h][y]
WHOMFLY[N_,h_][cd_CD]:=WHOMFLY[N][CD2CDP[cd]]*h^(Length[cd]/2)
WHOMFLY[N_][CDP[]]:=N
WHOMFLY[N_][cdp_CDP] := WHOMFLY[N][cdp] = Module[
  {m,chords,delta,delpol,CL},
  m=Length[cdp]/2;
  chords=CL @@ Union[Sort /@ Transpose[{Range[2m],List @@ cdp}]];
  delpol=Expand[Times @@ (chords /.
    {i_,j_} :> (delta[i-1,j]delta[j-1,i]-N*delta[i-1,i]delta[j-1,j]))
  ];
  delpol //. {
    delta[0,2m] -> N,
    delta[i_,i_] -> N,
    delta[i_,j_]delta[j_,k_] :> delta[i,k]
  }
]

Format[t[i_,j_],TeXForm] := StringForm["t^{`1``2`}",i,j]
Format[p[tl___],TeXForm] := StringJoin @@ ToString /@ TeXForm /@ {tl}
Format[b[x_,y_],TeXForm] := StringForm["[`1`,`2`]",
  ToString[TeXForm[x]],ToString[TeXForm[y]]]
Format[CD[cl___],TeXForm] :=
  StringJoin["[",StringJoin @@ ToString /@ {cl},"]"]

End[]
EndPackage[]
