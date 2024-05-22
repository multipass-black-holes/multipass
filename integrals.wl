(* ::Package:: *)

Slvc[x_,y_,z_]:=Piecewise[{{0,x<y},{1,x>y+z},{1/(Exp[z/(x-y)+z/(x-y-z)]+1),True}}]


LVCINT[y_?NumberQ] := NIntegrate[1/(1 + Exp[(1 - 2*x)/(x - x^2)]), {x, 0, y}]


Btilde[a_?NumberQ,b_?NumberQ,z_?NumberQ] := -Gamma[a] Gamma[b]/Gamma[a+b]+z^a Hypergeometric2F1[a,1-b,1+a,z]/a


mf1g[m_, a_, b_, d_, mgap_, mmin_, dm_]:=Piecewise[{
  {m^b*(1 + (2*a^2*m^(1/2)*(mgap - m)^(a - 1))/mgap^(a - 1/2))*Slvc[m, mmin, dm], m < mgap},
  {0, True}
}]


mf2g[m_, a_, b_, d_, mgap_, mmin_, dm_]:=Piecewise[{
  {0, m < mmin},
  {Slvc[m, mmin, dm], mmin < m < mgap + mmin + dm/2},
  {(m/(mgap + mmin + dm/2))^d, mgap + mmin + dm/2 < m}
}]


int0[mmin_, dm_, mgap_, a_, b_, d_, lam21_, lam12_, beta0_, beta1_]:=(
  (
    (1 + b - beta0)*mgap^(2 + 2*b)
    + (1 + b + beta0)*(dm + mmin)^(2 + 2*b)
    - 2*(1 + b)*mgap^(1 + b - beta0)*(dm + mmin)^(1 + b + beta0)
  )/(2*(1 + b)*(1 + b - beta0)*(1 + b + beta0))
  +Btilde[5/2 + 2*b, a, (dm + mmin)/mgap]
     *4*a^2*beta0*mgap^(2 + 2*b)/(1 + 2*b + b^2 - beta0^2)
  -Btilde[3/2 + b + beta0, a, (dm + mmin)/mgap]
     *2*a^2*mgap^(2 + 2*b)/(1 + b - beta0)
  +Btilde[3/2 + b - beta0, a, (dm + mmin)/mgap]
     *2*a^2*mgap^(1 + b - beta0)*(dm + mmin)^(1 + b + beta0)/(1 + b + beta0) 
  +Btilde[3/2 + b - beta0, a, (dm + mmin)/mgap]*Btilde[3/2 + b + beta0, a, (dm + mmin)/mgap]
     *4*a^4*mgap^(2 + 2*b)
  + 4*a^4*mgap^(1 + 2*b)*int[
       (1 - m1/mgap)^(-1 + a)*(m1/mgap)^(1/2 + b - beta0)*Btilde[3/2 + b + beta0, a, m1/mgap],
       {m1, dm + mmin, mgap}
    ]
  + (
     (mgap^(1 + b) - (dm + mmin)^(1 + b))/(1 + b)
     - 2*a^2*mgap^(1 + b)*Btilde[3/2 + b, a, (dm + mmin)/mgap]
   )*int[mf1g[m1, a, b, d, mgap, mmin, dm], {m1, mmin, dm + mmin}]
  + int[
    mf1g[m1, a, b, d, mgap, mmin, dm]*mf1g[m2, a, b, d, mgap, mmin, dm], 
    {m1, mmin, dm + mmin}, {m2, mmin, m1}
  ]
);


int1[mmin_, dm_, mgap_, a_, b_, d_, lam21_, lam12_, beta0_, beta1_]:=(
  4*a^2*Sqrt[mgap]*(dm + mmin)^(3/2 + b)*(1 - (dm + mmin)/mgap)^a/((3 + 2*a + 2*b)*(1 + beta1))
  + dm*(mgap^(1 + b) - (dm + mmin)^(1 + b))/(2*(1 + b))
  + (mgap^(2 + b) - (dm + mmin)^(2 + b))/((2 + b)*(1 + beta1))
  - (dm + mmin)^(1 + beta1)*(mgap^(1 + b - beta1) - (dm + mmin)^(1 + b - beta1))
     /((1 + b - beta1)*(1 + beta1))
  + a^2*mgap^(1 + b)*(-dm - (2*(3 + 2*b)*mgap)/((3 + 2*a + 2*b)*(1 + beta1)))
    *Btilde[3/2 + b, a, (dm + mmin)/mgap]
  + 2*a^2*mgap^(1 + b - beta1)*(dm + mmin)^(1 + beta1)/(1 + beta1)
    *Btilde[3/2 + b - beta1, a, (dm + mmin)/mgap]
  + dm*int[mf1g[m1, a, b, d, mgap, mmin, dm]*LVCINT[(m1 - mmin)/dm], {m1, mmin, dm + mmin}]
);


int2[mmin_, dm_, mgap_, a_, b_, d_, lam21_, lam12_, beta0_, beta1_]:=(
(4*a^2*Sqrt[mgap]*(dm + mmin)^(3/2 + b)*(-((dm - mgap + mmin)/mgap))^a)/((3 + 2*a + 2*b)*(-1 + beta1))
+ (
  (2*(mgap^(2 + b) - (dm + mmin)^(2 + b)))/(2 + b)
  + (2*(-(mgap^beta1*(dm + mmin)^(2 + b)) + mgap*(dm + mmin)^(1 + b + beta1)))/
    ((-1 + beta1)*mgap^beta1)
  + (mgap^(1 + b + beta1) - (dm + mmin)^(1 + b + beta1))*(
      2^beta1*d*(dm + 2*mmin)*(mgap*(dm + 2*mgap + 2*mmin))^beta1
      +2*mgap*(
        (-1 + beta1)*(dm + 2*mgap + 2*mmin)^(2*beta1)
        - d*(dm + 2*mgap + 2*mmin)^(2*beta1)
        + 2^beta1*d*(mgap*(dm + 2*mgap + 2*mmin))^beta1
      )
    )/((-1 + beta1)*(-1 + beta1 - d)*(mgap*(dm + 2*mgap + 2*mmin)^2)^beta1)
  )/(2*(1 + b + beta1))
- Btilde[3/2 + b, a, (dm + mmin)/mgap] * 2*a^2*(3 + 2*b)*mgap^(2 + b)/((3 + 2*a + 2*b)*(-1 + beta1))
- Btilde[3/2 + b + beta1, a, (dm + mmin)/mgap] *
    2^beta1*a^2*d*mgap^(1 + b + beta1)*(dm + 2*mgap + 2*mmin)^(1 - beta1)/(1 - 2*beta1 + beta1^2 + d - beta1*d)
+ int[mf1g[m1, a, b, d, mgap, mmin, dm], {m1, mmin, dm + mmin}](-2*dm - d*dm + 2*d*mgap - 2*mmin)/(2 + 2*d)
+ int[mf1g[m2, a, b, d, mgap, mmin, dm]*Slvc[m1, mmin, dm], {m1, mmin, dm + mmin}, {m2, mmin, m1}]
);


norms[mmin_?NumberQ, dm_?NumberQ, mgap_?NumberQ, a_?NumberQ, b_?NumberQ, d_?NumberQ, lam21_?NumberQ, lam12_?NumberQ, beta_?NumberQ]:={
  int0[mmin,dm,mgap,a,b,d,lam21,lam12,beta,beta],int1[mmin,dm,mgap,a,b,d,lam21,lam12,beta,beta],int2[mmin,dm,mgap,a,b,d,lam21,lam12,beta,beta]
}/.int->NIntegrate


(* ::Code:: *)
(*norms[3.98,5.62,70.2,0.169,-2.31,-5.58,1.85,-4.41,6.83]/{*)
(*  1.421332105126206`*^-3,*)
(*  2.382772068183056`*^-1,*)
(*  2.283400792009993`*)
(*}-1//Chop*)
