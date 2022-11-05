(* (c) 12/20/93, Dror Bar-Natan *)

{R[1]=1+t[1,2]/2, Phi[1]=1}

R[m_,0] := R[m-1]; Phi[m_,0] := Phi[m-1]

R[m_,1] := R[m,0]; Phi[m_,1] := Phi[m,1] =
  ModDegree[m+1,Phi[m,0]-Expand[(-1+NormSquared[Phi[m,0]])/2]]

RInverse[m_,i___] := RInverse[m,i] = ModDegree[m+1,Invert[R[m,i]]]
PhiInverse[m_,i___] := PhiInverse[m,i] = ModDegree[m+1,Invert[Phi[m,i]]]

psip[m_,i_] := ModDegree[m+1, -1 + Act[{1,2},3][RInverse[m,i]]**
  Phi[m,i]**Act[2,3][R[m,i]]**Act[1,3,2][PhiInverse[m,i]]**
  Act[1,3][R[m,i]]**Act[3,1,2][Phi[m,i]]]

psim[m_,i_] := ModDegree[m+1, -1 + Act[{1,2},3][R[m,i]]**Phi[m,i]**
  Act[2,3][RInverse[m,i]]**Act[1,3,2][PhiInverse[m,i]]**
  Act[1,3][RInverse[m,i]]**Act[3,1,2][Phi[m,i]]]

psi[m_,i_]:=Expand[(psip[m,i]+psim[m,i])/2]
psidiff[m_,i_]:=Expand[(psip[m,i]-psim[m,i])/2]

R[m_,2] := R[m,2] = ModDegree[m+1,
  R[m,1] - {APPower[t[1,2],m]}.
    (psidiff[m,1]~InTermsOf~{(d[1]~OX~1)[APPower[t[1,2],m]]})]
Phi[m_,2] := Phi[m,1]

R[m_,3] := R[m,2]; Phi[m_,3] := Phi[m,2]

R[m_,4] := R[m,3]; Phi[m_,4] := Phi[m,4] = Phi[m,3]-Expand[psi[m,3]/3]

mu[m_,i_] := mu[m,i] = ModDegree[m+1, -1 + Phi[m,i]**
  Act[1,{2,3},4][Phi[m,i]]**Act[2,3,4][Phi[m,i]]**
  Act[1,2,{3,4}][PhiInverse[m,i]]**Act[{1,2},3,4][PhiInverse[m,i]]]

R[m_,5] := R[m,4]; Phi[m_,5] := Phi[m,5] = Phi[m,4] + Expand[
  HarrisonBasis[3,m].((-mu[m,4])~InTermsOf~(d[3] /@ HarrisonBasis[3,m]))]

R[m_] := R[m,5]; Phi[m_] := Phi[m,5]

ZInfinity[m_] := ZInfinity[m] = AP2CD[{3,2,1},Phi[m]]
ZInfinityInverse[m_]:=ZInfinityInverse[m]=ModDegree[m+1,Invert[ZInfinity[m]]]

Trefoil1[m_] := Trefoil1[m]=ModDegree[m+1,AP2CD[{1,3,4,2},
  Act[1,2,{3,4}][Phi[m]]**
  Act[2,3,4][PhiInverse[m]]**
  Act[2,3][RInverse[m]]**
  Act[2,3][RInverse[m]]**
  Act[2,3][RInverse[m]]**
  Act[3,2,4][Phi[m]]**
  Act[1,3,{2,4}][PhiInverse[m]]
]**ZInfinityInverse[m]**ZInfinityInverse[m]]

Trefoil2[m_] := Trefoil2[m]=ModDegree[m+1,AP2CD[{1,3,2,4},
  PhiInverse[m]**
  Act[{1,2},3,4][Phi[m]]**
  RInverse[m]**
  Act[3,4][RInverse[m]]**
  Act[2,1,{3,4}][Phi[m]]**
  Act[1,4,3][PhiInverse[m]]**
  Act[1,4][R[m]]**
  Act[4,1,3][Phi[m]]
]**ZInfinityInverse[m]**ZInfinityInverse[m]]
