(* (c) 12/20/93, Dror Bar-Natan *)

BeginPackage["LinAlg`"]

SpanAppend::usage = "SpanAppend[plist,p] returns a list of polynomials \
	whose span is the same as the span of Append[plist,p]. plist \
	should have itself been made by SpanAppend."

SpanningIndices::usage = "SpanningIndices[plist] returns a list of indices \
	of elements in plist that linearly span plist"

LinearSpan::usage = "LinearSpan[plist] returns a basis to the linear span of \
	the polynomials in plist."

LinearDependencies::usage = "LinearDependencies[plist] returns a basis to \
	the space of linear dependencies between the polynomials in plist."

InTermsOf::usage = "p~InTermsOf~plist returns a list of numbers l for which \
	p is equal to l.plist, if such a list exists. Otherwise it returns \
	junk."
InTermsOf::NotInSpan = "Vector `1` is not in the linear span of `2`"

Begin["`private`"]

SpanningIndices[plist_List] := Red[plist, ReturnValue -> BasisAt]

LinearSpan[plist_List] := Red[plist, ReturnValue -> UnreducedBasis]

LinearDependencies[plist_List] := Red[plist, ReturnValue -> Kernel]

InTermsOf[vec_,bas_]:= Module[{mat,at},
  {mat,at}=Red[Append[bas,vec],ReturnValue -> {ReductionMatrix,BasisAt}];
  If[Last[at] > Length[bas],Message[InTermsOf::NotInSpan,vec,bas]];
  -Drop[Last[mat],-1]
]
	
LastTerm[exp_Plus]:=Last[exp];
LastTerm[exp_]:=exp;

FirstTerm[exp_Plus]:=First[exp];
FirstTerm[exp_]:=exp;

Options[Red] := {ReturnValue -> ReductionMatrix, Expander -> Expand,
		 Printing -> False}
Red[l_,opts___]:=Module[{v={},i,t,j,at={},ll,mat,coeff,lt,
    retval = (ReturnValue /. {opts} /. Options[Red]),
    expand = (Expander /. {opts} /. Options[Red]),
    printing = (Printing /. {opts} /. Options[Red])
  },
  ll=Length[l];
  If[printing,Print["Length[l]      ",ll]];
  mat=IdentityMatrix[ll];
  For[i=1,i<=ll,++i,(
    t=expand[l[[i]]];
    For[j=1,j<=Length[v],++j,(
      coeff = -Coefficient[t,LastTerm[v[[j]]]];
      If[!(coeff===0),
        mat[[i]] += coeff*mat[[at[[j]]]];
        t=Expand[t+v[[j]]*coeff];
      ]
    )];
    If[!(t===0),(
      lt=LastTerm[Expand[t]];
      coeff=If[Head[lt]===Times,Times @@ Select[List @@ lt,NumberQ],1];
      mat[[i]] /= coeff;
      t=Expand[t/coeff];
      AppendTo[v,t];
      AppendTo[at,i];
      If[printing,Print[i," appended"]];
    ),If[printing,Print[i," discarded"]]
    ];
  )];
  If[printing,Print["Length[v]      ",Length[v]]];
  If[printing,Print["at             ",at]];
  retval /. {
    ReductionMatrix :> mat,
    BasisAt :> at,
    Basis :> v,
    UnreducedBasis :> expand /@ l[[at]],
    Kernel :> mat[[Complement[Range[ll],at]]]
  }
]

red[l_,opts___]:=Module[{v,lv=0,i,t,j,at={},ll,mat,coeff,lt,
    retval = (ReturnValue /. {opts} /. Options[Red]),
    expand = (Expander /. {opts} /. Options[Red]),
    printing = (Printing /. {opts} /. Options[Red])
  },
  ll=Length[l];
  If[printing,Print["Length[l]      ",ll]];
  mat=IdentityMatrix[ll];
  For[i=1,i<=ll,++i,(
    t=expand[l[[i]]];
    For[j=1,j<=lv,++j,(
      coeff = -Coefficient[t,LastTerm[v[j]]];
      If[!(coeff===0),
        mat[[i]] += coeff*mat[[at[[j]]]];
        t=Expand[t+v[j]*coeff];
      ]
    )];
    If[!(t===0),(
      lt=LastTerm[Expand[t]];
      coeff=If[Head[lt]===Times,Times @@ Select[List @@ lt,NumberQ],1];
      mat[[i]] /= coeff;
      t=Expand[t/coeff];
      ++lv;
      v[lv]=t;
      AppendTo[at,i];
      If[printing,Print[i," appended"]];
    ),If[printing,Print[i," discarded"]]
    ];
  )];
  If[printing,Print["Length[v]      ",lv]];
  If[printing,Print["at             ",at]];
  retval /. {
    Basis :> v,
    ReductionMatrix :> mat,
    BasisAt :> at,
    Kernel :> mat[[Complement[Range[ll],at]]]
  }
]

Clear[NZLength,ZDrop]
NZLength[l_List]:=Length[Select[l,(!(#===0))&]]
ZDrop[l_List]:=Select[l,(!(#===0))&]

NumericalCoefficient[c_?NumberQ]:=c
NumericalCoefficient[c_?NumberQ * _] := c
NumericalCoefficient[_]=1

Normalize[expr_] := Expand[expr/NumericalCoefficient[FirstTerm[expr]]]

SpanAppend[l_List,t_]:=Module[{a=t,b,r},
  If[a=!=0,a=Normalize[a]];
  Join[
    Function[b,
      If[a===0,b,
      If[FirstTerm[a]=!=(ft=FirstTerm[b]),
        a=Expand[a-Coefficient[a,ft]*b];
        b,
        If[If[Head[a]===Plus,Length[a],1] < If[Head[b]===Plus,Length[b],1],
          r=a;
          a=Expand[b-Coefficient[b,FirstTerm[a]]*a]
        ,
          r=b;
          a=Expand[a-Coefficient[a,ft]*b]
        ];
        If[a=!=0,a=Normalize[a]];
        r
      ]]
    ] /@ l,
    If[a===0,{},{a}]
  ]
]

End[]
EndPackage[]
