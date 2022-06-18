(* ::Package:: *)

(* ::Title:: *)
(*grid and integration *)


(* ::Section:: *)
(*1. Global constants and basic functions*)


(* ::Input::Initialization:: *)
hbar=1.0545718176461565*10^-27;
e=4.803204712570263*10^(-10);
me=9.1093837015*10^(-28);
mp=1.67262192369*10^-24;
A=1;Z=1;Ye=1;
kb=1.380649*10^(-16);
kev2erg=1.6021766339999998*10^-9;
c=29979245800.0;
\[Sigma]T=8*\[Pi]/3*e^4/(c^4 me^2);
\[Kappa]T=\[Sigma]T/mp;
G=6.674299999999999*10^-8;
Ms=1.988409870698051*10^33;
km=10^5;
\[Alpha]F=e^2/(hbar c);
BQ=me^2*c^3/e/hbar;
logspace[a_,b_,n_]:=10.0^Range[a,b,(b-a)/(n-1)];
linearspace[a_,b_,n_]:=Range[a,b,(b-a)/(n-1)];
C0[a_]:=1/a-Exp[a]*ExpIntegralE[1,a];
C1[a_]:=(1+a)*Exp[a]*ExpIntegralE[1,a]-1;
radian[x_]:=x/180*\[Pi];
Planck[Es_,T_]:=2*hbar*2*\[Pi]*(Es*kev2erg/hbar/2/\[Pi])^3/c^2*1/(Exp[(Es*kev2erg)/(kb*T)]-1);


(* ::Section:: *)
(*Temperature profile*)


(* ::Input::Initialization:: *)
temperature[\[Tau]_]:=Module[{\[Tau]mid,x},
\[Tau]mid=4.27;
a1=-0.0599;
a2=0.192;
a3=0.0225;
a4=0.0115;
a5= -0.0072;
a6=0.00116;
b3=0.109;
b4=0.0828;
b5=0.0256;
b6=0.00286;
x=Log[10,\[Tau]]-Log[10,\[Tau]mid];
Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]>\[Tau]mid && \[Tau]<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]<= \[Tau]mid && \[Tau]>= 10^-3}}]];
LogLogPlot[temperature[\[Tau]],{\[Tau],0.001,10^4},Background->White,PlotTheme->"Monochrome"]


(* ::Section:: *)
(*Scatter opacity*)


(* ::Input::Initialization:: *)
Aint[Es_,B_,\[Tau]_]:=Module[{Ec,\[Omega],\[Omega]c,Ece,Eci,Epe,Epi,ue,ui,ve,vi,b,a,q,m,ri,\[Epsilon]p,\[Eta]p,g,\[Epsilon],\[Eta],\[Beta]i,K1i,K2i,Kz1i,Kz2i,ep1i,em1i,ep2i,em2i,e01i,e02i,Ap1,Am1,A01,Ap2,Am2,A02,Ap,Am,A0,\[Rho],x,T,Mns,Rns,gns,\[Tau]mid,\[Delta]v},
\[Tau]mid=4.27;
a1=-0.0599;
a2=0.192;
a3=0.0225;
a4=0.0115;
a5= -0.0072;
a6=0.00116;
b3=0.109;
b4=0.0828;
b5=0.0256;
b6=0.00286;
x=Log[10,\[Tau]]-Log[10,\[Tau]mid];
T=Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]>\[Tau]mid && \[Tau]<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]<= \[Tau]mid && \[Tau]>= 10^-3}}];
Mns=1.4*Ms;
Rns=10^6;
gns=(G*Mns)/Rns^2 (1-(2G Mns)/(Rns c^2))^(-1/2);
\[Rho]=mp*gns/(2*\[Kappa]T*kb*T)*\[Tau];
Ec=Es*kev2erg;
\[Omega]=Ec/hbar;
\[Omega]c=(e B)/(me c);
Ece=hbar*\[Omega]c/kev2erg;Eci=Ece*(Z me)/(A mp);
Epe=hbar*(4*\[Pi]*\[Rho]/mp * Ye * e^2/me)^0.5/kev2erg;
Epi=((Z me)/(A mp))^0.5*Epe;
ue=(Ece/Es)^2;ui=(Eci/Es)^2;ve=(Epe/Es)^2;vi=(Epi/Es)^2;

b=B/BQ;
\[Delta]v=\[Alpha]F/(45 \[Pi]) b^2;
(*turn on the vacuum*)
a=1-2*\[Delta]v;
q=7\[Delta]v;
m=-4\[Delta]v;
(*turn off the vacuum*)
(*a=1;
q=0;
m=0;*)
\[Epsilon]p=1-(ve(1-(A mp)/(Z me) ui))/((1-ui)(1-ue));
\[Eta]p=1-ve;
g=(ve*ue^0.5)/((1-ui)(1-ue));
\[Epsilon]=\[Epsilon]p+a-1;
\[Eta]=\[Eta]p+a+q-1;
ri=1+m/a Sin[\[Theta]i]^2;
\[Beta]i=-((\[Epsilon]^2-g^2-\[Epsilon] *\[Eta](1+m/a))/(2g \[Eta]))*Sin[\[Theta]i]^2/Cos[\[Theta]i];
K1i=\[Beta]i(1-(1+ri/\[Beta]i^2)^(1/2));
K2i=\[Beta]i(1+(1+ri/\[Beta]i^2)^(1/2));
Kz1i=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]i]Cos[\[Theta]i]K1i + g * Sin[\[Theta]i])/(\[Epsilon]p*Sin[\[Theta]i]^2+(\[Eta]p+q)Cos[\[Theta]i]^2+a-1);

Kz2i=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]i]Cos[\[Theta]i]K2i + g * Sin[\[Theta]i])/(\[Epsilon]p*Sin[\[Theta]i]^2+(\[Eta]p+q)Cos[\[Theta]i]^2+a-1);
ep1i=(1+K1i*Cos[\[Theta]i]+Kz1i*Sin[\[Theta]i])^2/(2(1+K1i^2+Kz1i^2));
ep2i=(1+K2i*Cos[\[Theta]i]+Kz2i*Sin[\[Theta]i])^2/(2(1+K2i^2+Kz2i^2));
em1i=(1-(K1i*Cos[\[Theta]i]+Kz1i*Sin[\[Theta]i]))^2/(2(1+K1i^2+Kz1i^2));
em2i=(1-(K2i*Cos[\[Theta]i]+Kz2i*Sin[\[Theta]i]))^2/(2(1+K2i^2+Kz2i^2));
e01i=(K1i * Sin[\[Theta]i]-Kz1i*Cos[\[Theta]i])^2/(1+K1i^2+Kz1i^2);
e02i=(K2i * Sin[\[Theta]i]-Kz2i*Cos[\[Theta]i])^2/(1+K2i^2+Kz2i^2);
Ap=NIntegrate[3/4*Sin[\[Theta]i]*(ep1i+ep2i),{\[Theta]i,0,\[Pi]},WorkingPrecision->10,PrecisionGoal->6];
Am=NIntegrate[3/4*Sin[\[Theta]i]*(em1i+em2i),{\[Theta]i,0,\[Pi]},WorkingPrecision->10,PrecisionGoal->6];
A0=NIntegrate[3/4*Sin[\[Theta]i]*(e01i+e02i),{\[Theta]i,0,\[Pi]},WorkingPrecision->10,PrecisionGoal->6];
Ap2=NIntegrate[3/4*Sin[\[Theta]i]*ep2i,{\[Theta]i,0,\[Pi]},WorkingPrecision->10,PrecisionGoal->6];Am2=NIntegrate[3/4*Sin[\[Theta]i]*em2i,{\[Theta]i,0,\[Pi]},WorkingPrecision->10,PrecisionGoal->6];
A02=NIntegrate[3/4*Sin[\[Theta]i]*e02i,{\[Theta]i,0,\[Pi]},WorkingPrecision->10,PrecisionGoal->6];
Ap1=Ap-Ap2;
Am1=Am-Am2;
A01=A0-A02;
{Ap1,Ap2,Am1,Am2,A01,A02}
]


(* ::Input::Initialization:: *)
scatteropacity[Es_,B_,\[Tau]_,\[Theta]b_,Ap1_,Ap2_,Am1_,Am2_,A01_,A02_]:=Module[{Ec,\[Omega],\[Omega]c,Ece,Eci,Epe,Epi,ue,ui,ve,vi,b,a,q,m,r,\[Epsilon]p,\[Eta]p,g,\[Epsilon],\[Eta],\[Beta],K1,K2,Kz1,Kz2,ep1,em1,ep2,em2,e01,e02,ne,ni,gl,gp,y1,y2,\[Omega]ci,\[Alpha]0,\[Nu]re,\[Nu]ep,\[Nu]em,\[Nu]e0,\[Nu]ip,\[Nu]im,\[Nu]i0,sc11,sc12,sc21,sc22,\[Tau]mid,x,T,\[Rho],Mns,Rns,gns,\[Delta]v},
\[Tau]mid=4.27;
a1=-0.0599;
a2=0.192;
a3=0.0225;
a4=0.0115;
a5= -0.0072;
a6=0.00116;
b3=0.109;
b4=0.0828;
b5=0.0256;
b6=0.00286;
x=Log[10,\[Tau]]-Log[10,\[Tau]mid];
T=Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]>\[Tau]mid && \[Tau]<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]<= \[Tau]mid && \[Tau]>= 10^-3}}];
Mns=1.4*Ms;
Rns=10^6;
gns=(G*Mns)/Rns^2 (1-(2G Mns)/(Rns c^2))^(-1/2);
\[Rho]=mp*gns/(2*\[Kappa]T*kb*T)*\[Tau];
Ec=Es*kev2erg;
\[Omega]=Ec/hbar;
\[Omega]c=(e B)/(me c);
\[Omega]ci= (Z e B)/(A mp c);
Ece=hbar*\[Omega]c/kev2erg;Eci=Ece*(Z me)/(A mp);
Epe=hbar*(4*\[Pi]*\[Rho]/mp * Ye * e^2/me)^0.5/kev2erg;
Epi=((Z me)/(A mp))^0.5*Epe;
ue=(Ece/Es)^2;ui=(Eci/Es)^2;ve=(Epe/Es)^2;vi=(Epi/Es)^2;
b=B/BQ;
(*turn on the vacuum*)
\[Delta]v=\[Alpha]F/(45 \[Pi]) b^2;
a=1-2*\[Delta]v;
q=7\[Delta]v;
m=-4\[Delta]v;
(*turn off the vacuum*)
(*a=1;
q=0;
m=0;*)
r=1+m/a Sin[\[Theta]b]^2;
\[Epsilon]p=1-(ve(1-(A mp)/(Z me) ui))/((1-ui)(1-ue));
\[Eta]p=1-ve;
g=(ve*ue^0.5)/((1-ui)(1-ue));
\[Epsilon]=\[Epsilon]p+a-1;
\[Eta]=\[Eta]p+a+q-1;
\[Beta]=-((\[Epsilon]^2-g^2-\[Epsilon] *\[Eta](1+m/a))/(2g \[Eta]))*Sin[\[Theta]b]^2/Cos[\[Theta]b];
K1=\[Beta](1-(1+r/\[Beta]^2)^(1/2));
K2=\[Beta](1+(1+r/\[Beta]^2)^(1/2));
Kz1=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K1 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1);

Kz2=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K2 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1);
ep1=(1+K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b])^2/(2(1+K1^2+Kz1^2));
ep2=(1+K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b])^2/(2(1+K2^2+Kz2^2));
em1=(1-(K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b]))^2/(2(1+K1^2+Kz1^2));
em2=(1-(K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b]))^2/(2(1+K2^2+Kz2^2));
e01=(K1 * Sin[\[Theta]b]-Kz1*Cos[\[Theta]b])^2/(1+K1^2+Kz1^2);
e02=(K2 * Sin[\[Theta]b]-Kz2*Cos[\[Theta]b])^2/(1+K2^2+Kz2^2);
ne=\[Rho]/mp;
ni=ne;
y1=hbar*\[Omega]c/(kb*T);
y2=\[Omega]/\[Omega]c;
gl=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*C1[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
gp=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*2*y2*Exp[2x1]C0[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
\[Alpha]0=4*\[Pi]^2 Z^2 \[Alpha]F^3 (hbar^2 c^2)/me^2 (2 me/(\[Pi] kb T))^(1/2) (ne ni)/\[Omega]^3 (1-Exp[(-hbar \[Omega])/(kb *T)]);
\[Nu]re=(2*e^2)/(3 me c^3) \[Omega]^2;
\[Nu]ep=\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re;
\[Nu]em=\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re;
\[Nu]e0=\[Nu]re+(\[Alpha]0*gp)/(ne \[Sigma]T) \[Nu]re;
\[Nu]ip=(\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re) me/mp;
\[Nu]im=(\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re) me/mp;
\[Nu]i0=(\[Nu]re+(\[Alpha]0*gp)/(ne \[Sigma]T) \[Nu]re) me/mp;

sc11=(ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]+1*\[Omega]c)^2+\[Nu]ep^2))ep1*Ap1+(\[Omega]^2/((\[Omega]-1*\[Omega]c)^2+\[Nu]em^2))em1*Am1+(\[Omega]^2/((\[Omega]+0*\[Omega]c)^2+\[Nu]e0^2))e01*A01)+(me/ mp)^2 (ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]-1*\[Omega]ci)^2+\[Nu]ip^2))ep1*Ap1+(\[Omega]^2/((\[Omega]+1*\[Omega]ci)^2+\[Nu]im^2))em1*Am1+(\[Omega]^2/((\[Omega]+0*\[Omega]ci)^2+\[Nu]i0^2))e01*A01);

sc12=(ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]+1*\[Omega]c)^2+\[Nu]ep^2))ep1*Ap2+(\[Omega]^2/((\[Omega]-1*\[Omega]c)^2+\[Nu]em^2))em1*Am2+(\[Omega]^2/((\[Omega]+0*\[Omega]c)^2+\[Nu]e0^2))e01*A02)+(me/ mp)^2 (ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]-1*\[Omega]ci)^2+\[Nu]ip^2))ep1*Ap2+(\[Omega]^2/((\[Omega]+1*\[Omega]ci)^2+\[Nu]im^2))em1*Am2+(\[Omega]^2/((\[Omega]+0*\[Omega]ci)^2+\[Nu]i0^2))e01*A02);
sc21=(ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]+1*\[Omega]c)^2+\[Nu]ep^2))ep2*Ap1+(\[Omega]^2/((\[Omega]-1*\[Omega]c)^2+\[Nu]em^2))em2*Am1+(\[Omega]^2/((\[Omega]+0*\[Omega]c)^2+\[Nu]e0^2))e02*A01)+(me/ mp)^2 (ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]-1*\[Omega]ci)^2+\[Nu]ip^2))ep2*Ap1+(\[Omega]^2/((\[Omega]+1*\[Omega]ci)^2+\[Nu]im^2))em2*Am1+(\[Omega]^2/((\[Omega]+0*\[Omega]ci)^2+\[Nu]i0^2))e02*A01);

sc22=(ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]+1*\[Omega]c)^2+\[Nu]ep^2))ep2*Ap2+(\[Omega]^2/((\[Omega]-1*\[Omega]c)^2+\[Nu]em^2))em2*Am2+(\[Omega]^2/((\[Omega]+0*\[Omega]c)^2+\[Nu]e0^2))e02*A02)+(me/ mp)^2 (ne \[Sigma]T)/\[Rho] ((\[Omega]^2/((\[Omega]-1*\[Omega]ci)^2+\[Nu]ip^2))ep2*Ap2+(\[Omega]^2/((\[Omega]+1*\[Omega]ci)^2+\[Nu]im^2))em2*Am2+(\[Omega]^2/((\[Omega]+0*\[Omega]ci)^2+\[Nu]i0^2))e02*A02);
{sc11,sc12,sc21,sc22}

]


(* ::Section:: *)
(*free free opacity*)


(* ::Input::Initialization:: *)
freeopacity[Es_,B_,\[Tau]_,\[Theta]b_]:=Module[{Ec,\[Omega],\[Omega]c,Ece,Eci,Epe,Epi,ue,ui,ve,vi,b,a,q,m,r,\[Epsilon]p,\[Eta]p,g,\[Epsilon],\[Eta],\[Beta],K1,K2,Kz1,Kz2,ep1,em1,ep2,em2,e01,e02,ne,ni,\[Alpha]0,gl,gp,y1,y2,\[Omega]ci,\[Nu]re,\[Nu]ep,\[Nu]em,\[Nu]e0,\[Nu]ip,\[Nu]im,\[Nu]i0,ffe1,ffi1,ffe2,ffi2,ff1,ff2,T,\[Rho],Mns,Rns,gns,\[Tau]mid,x,\[Delta]v},
\[Tau]mid=4.27;
a1=-0.0599;
a2=0.192;
a3=0.0225;
a4=0.0115;
a5= -0.0072;
a6=0.00116;
b3=0.109;
b4=0.0828;
b5=0.0256;
b6=0.00286;
x=Log[10,\[Tau]]-Log[10,\[Tau]mid];
T=Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]>\[Tau]mid && \[Tau]<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]<= \[Tau]mid && \[Tau]>= 10^-3}}];
Mns=1.4*Ms;
Rns=10^6;
gns=(G*Mns)/Rns^2 (1-(2G Mns)/(Rns c^2))^(-1/2);
\[Rho]=mp*gns/(2*\[Kappa]T*kb*T)*\[Tau];
Ec=Es*kev2erg;
\[Omega]=Ec/hbar;
\[Omega]c=(e B)/(me c);
\[Omega]ci= (Z e B)/(A mp c);
Ece=hbar*\[Omega]c/kev2erg;Eci=Ece*(Z me)/(A mp);
Epe=hbar*(4*\[Pi]*\[Rho]/mp * Ye * e^2/me)^0.5/kev2erg;
Epi=((Z me)/(A mp))^0.5*Epe;
ue=(Ece/Es)^2;ui=(Eci/Es)^2;ve=(Epe/Es)^2;vi=(Epi/Es)^2;
b=B/BQ;
(*turn on the vacuum*)
\[Delta]v=\[Alpha]F/(45 \[Pi]) b^2;
a=1-2*\[Delta]v;
q=7\[Delta]v;
m=-4\[Delta]v;
(*turn off the vacuum*)
(*a=1;
q=0;
m=0;*)
r=1+m/a Sin[\[Theta]b]^2;
\[Epsilon]p=1-(ve(1-(A mp)/(Z me) ui))/((1-ui)(1-ue));
\[Eta]p=1-ve;
g=(ve*ue^0.5)/((1-ui)(1-ue));
\[Epsilon]=\[Epsilon]p+a-1;
\[Eta]=\[Eta]p+a+q-1;
\[Beta]=-((\[Epsilon]^2-g^2-\[Epsilon] *\[Eta](1+m/a))/(2g \[Eta]))*Sin[\[Theta]b]^2/Cos[\[Theta]b];
K1=\[Beta](1-(1+r/\[Beta]^2)^(1/2));
K2=\[Beta](1+(1+r/\[Beta]^2)^(1/2));
Kz1=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K1 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1);

Kz2=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K2 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1);
ep1=(1+K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b])^2/(2(1+K1^2+Kz1^2));
ep2=(1+K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b])^2/(2(1+K2^2+Kz2^2));
em1=(1-(K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b]))^2/(2(1+K1^2+Kz1^2));
em2=(1-(K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b]))^2/(2(1+K2^2+Kz2^2));
e01=(K1 * Sin[\[Theta]b]-Kz1*Cos[\[Theta]b])^2/(1+K1^2+Kz1^2);
e02=(K2 * Sin[\[Theta]b]-Kz2*Cos[\[Theta]b])^2/(1+K2^2+Kz2^2);
ne=\[Rho]/mp;
ni=ne;
\[Alpha]0=4*\[Pi]^2 Z^2 \[Alpha]F^3 (hbar^2 c^2)/me^2 (2 me/(\[Pi] kb T))^(1/2) (ne ni)/\[Omega]^3 (1-Exp[(-hbar \[Omega])/(kb *T)]);
y1=hbar*\[Omega]c/(kb*T);
y2=\[Omega]/\[Omega]c;
gl=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*C1[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
gp=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*2*y2*Exp[2x1]C0[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
\[Nu]re=(2*e^2)/(3 me c^3) \[Omega]^2;
\[Nu]ep=\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re;
\[Nu]em=\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re;
\[Nu]e0=\[Nu]re+(\[Alpha]0*gp)/(ne \[Sigma]T) \[Nu]re;
\[Nu]ip=(\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re) me/mp;
\[Nu]im=(\[Nu]re+(\[Alpha]0*gl)/(ne \[Sigma]T) \[Nu]re) me/mp;
\[Nu]i0=(\[Nu]re+(\[Alpha]0*gp)/(ne \[Sigma]T) \[Nu]re) me/mp;
ffe1=\[Alpha]0/\[Rho] ((\[Omega]^2/((\[Omega]+1*\[Omega]c)^2+\[Nu]ep^2))ep1*gl+(\[Omega]^2/((\[Omega]-1*\[Omega]c)^2+\[Nu]em^2))em1*gl+(\[Omega]^2/((\[Omega]+0*\[Omega]c)^2+\[Nu]e0^2))e01*gp);
ffi1=1/Z^3 ((Z^2 me)/(A mp))^2 \[Alpha]0/\[Rho] ((\[Omega]^2/((\[Omega]-1*\[Omega]ci)^2+\[Nu]ip^2))ep1*gl+(\[Omega]^2/((\[Omega]+1*\[Omega]ci)^2+\[Nu]im^2))em1*gl+(\[Omega]^2/((\[Omega]+0*\[Omega]ci)^2+\[Nu]i0^2))e01*gp);
ffe2=\[Alpha]0/\[Rho] ((\[Omega]^2/((\[Omega]+1*\[Omega]c)^2+\[Nu]ep^2))ep2*gl+(\[Omega]^2/((\[Omega]-1*\[Omega]c)^2+\[Nu]em^2))em2*gl+(\[Omega]^2/((\[Omega]+0*\[Omega]c)^2+\[Nu]e0^2))e02*gp);
ffi2=1/Z^3 ((Z^2 me)/(A mp))^2 \[Alpha]0/\[Rho] ((\[Omega]^2/((\[Omega]-1*\[Omega]ci)^2+\[Nu]ip^2))ep2*gl+(\[Omega]^2/((\[Omega]+1*\[Omega]ci)^2+\[Nu]im^2))em2*gl+(\[Omega]^2/((\[Omega]+0*\[Omega]ci)^2+\[Nu]i0^2))e02*gp);


ff1=ffe1+ffi1;
ff2=ffe2+ffi2;
{ff1,ff2}

]


(* ::Section:: *)
(*jump condition*)


(* ::Input::Initialization:: *)
\[Tau]reson[Es_,B_]:=Module[{q,m,\[Omega]c,Ebi,Mns,Rns,Ead,b,\[Delta]v,fB,gns,\[Rho]v,\[Tau]v},
b=B/BQ;
\[Delta]v=\[Alpha]F/(45*\[Pi]) b^2;
q=7\[Delta]v;
m=-4\[Delta]v;
fB=((3*\[Delta]v)/(q+m))^0.5;
Mns=1.4*Ms;
Rns=10^6;
gns=(G*Mns)/Rns^2 (1-(2G Mns)/(Rns c^2))^(-1/2);
\[Rho]v=0.96* Ye^-1 (Es/1)^2 (B/10^14)^2 fB^-2;
\[Tau]v=\[Tau]x/.FindRoot[temperature[\[Tau]x]-\[Tau]x*(mp gns)/(\[Kappa]T 2 \[Rho]v kb)==0 ,{\[Tau]x,((mp gns)/(\[Kappa]T 2 \[Rho]v kb*5*10^6))^-1}];
\[Tau]v]

Pjump[Es_,B_,\[Theta]b_]:=Module[{\[Tau]v,Tv,H\[Rho],jump,b,q,m,\[Delta]v,fB,Mns,Rns,gns,\[Omega]c,Ebi,\[Rho]v,Ead},
b=B/BQ;
\[Delta]v=\[Alpha]F/(45*\[Pi]) b^2;
q=7\[Delta]v;
m=-4\[Delta]v;
fB=((3*\[Delta]v)/(q+m))^0.5;
Mns=1.4*Ms;
Rns=10^6;
gns=(G*Mns)/Rns^2 (1-(2G Mns)/(Rns c^2))^(-1/2);
\[Omega]c=(e B)/(me c);
Ebi=hbar*\[Omega]c*(me/mp)/kev2erg;
\[Rho]v=0.96* Ye^-1 (Es/1)^2 (B/10^14)^2 fB^-2;
\[Tau]v=\[Tau]x/.FindRoot[temperature[\[Tau]x]-\[Tau]x*(mp gns)/(\[Kappa]T 2 \[Rho]v kb)==0 ,{\[Tau]x,((mp gns)/(\[Kappa]T 2 \[Rho]v kb*5*10^6))^-1}];
Tv=temperature[\[Tau]v];

H\[Rho]=Abs[2*kb*Tv/(mp * gns *Cos[\[Theta]b])];
Ead=2.52*((fB*Tan[\[Theta]b])^2)^(1/3)*(Abs[1-(Ebi/Es)^2])^(2/3) (1/H\[Rho])^(1/3);
jump=Exp[-\[Pi]/2 (Es/Ead)^3];
N[jump]
]


(* ::Section:: *)
(*Generating the grid*)


(* ::Input::Initialization:: *)
Es=3;
B=4*10^13;
\[Tau]v=\[Tau]reson[Es,B];
\[Tau]x1=logspace[-3,Log[10,\[Tau]v],(IntegerPart[Log[10,\[Tau]v]]+3)*20];
\[Tau]x2=logspace[Log[10,\[Tau]v],Log[10,2]+4,IntegerPart[(Log[10,2]+4-Log[10,\[Tau]v])]*60];
\[Tau]=Reverse[DeleteDuplicates[Flatten[{\[Tau]x1,\[Tau]x2}]]];
position=IntegerPart[(Log[10,2]+4-Log[10,\[Tau]v])]*60;
\[Theta]x=N[Range[0.001,89.9,(89.9-0.001)/150]]/180*\[Pi];
I1u=Table[0,Length[\[Tau]],Length[\[Theta]x]];
I2u=Table[0,Length[\[Tau]],Length[\[Theta]x]];
I1d=Table[0,Length[\[Tau]],Length[\[Theta]x]];
I2d=Table[0,Length[\[Tau]],Length[\[Theta]x]];
I1u[[1]]=Flatten[Table[1/2,1,Length[\[Theta]x]]*Planck[Es,temperature[First[\[Tau]]]]];
I2u[[1]]=Flatten[Table[1/2,1,Length[\[Theta]x]]*Planck[Es,temperature[First[\[Tau]]]]];
I1d[[Length[\[Tau]]]]=Flatten[Table[0,1,Length[\[Theta]x]]];
I2d[[Length[\[Tau]]]]=Flatten[Table[0,1,Length[\[Theta]x]]];
(*give the initial data for source function*)
S1=Table[0,Length[\[Tau]],Length[\[Theta]x]];
S2=Table[0,Length[\[Tau]],Length[\[Theta]x]];

u1=Table[0,Length[\[Tau]]];
u2=Table[0,Length[\[Tau]]];
\[Delta]\[Tau]=Table[0,Length[\[Tau]]];
\[Delta]\[Tau]1=Table[0,Length[\[Tau]]];
For[j=1,j<=Length[\[Theta]x],j++,
For[
i=1,i<=Length[\[Tau]],i++,S1[[i,j]]=1/2*Planck[Es,temperature[\[Tau][[i]]]];S2[[i,j]]=1/2*Planck[Es,temperature[\[Tau][[i]]]]];
]


(* ::Input::Initialization:: *)
ff1=Table[0,{Length[\[Tau]]},{Length[\[Theta]x]}];
ff2=Table[0,{Length[\[Tau]]},{Length[\[Theta]x]}];
sc11=Table[0,{Length[\[Tau]]},{Length[\[Theta]x]}];
sc12=Table[0,{Length[\[Tau]]},{Length[\[Theta]x]}];
sc21=Table[0,{Length[\[Tau]]},{Length[\[Theta]x]}];
sc22=Table[0,{Length[\[Tau]]},{Length[\[Theta]x]}];
\[Kappa]t1=Table[0,Length[\[Tau]],Length[\[Theta]x]];
\[Kappa]t2=Table[0,Length[\[Tau]],Length[\[Theta]x]];
For[i=1,i<=Length[\[Tau]],i++,
{Ap1,Ap2,Am1,Am2,A01,A02}=Aint[Es,B,\[Tau][[i]]];
For[j=1,j<=Length[\[Theta]x],j++,
free=freeopacity[Es,B,\[Tau][[i]],\[Theta]x[[j]]];
ff1[[i,j]]=free[[1]];ff2[[i,j]]=free[[2]];
scatter=scatteropacity[Es,B,\[Tau][[i]],\[Theta]x[[j]],Ap1,Ap2,Am1,Am2,A01,A02];
sc11[[i,j]]=scatter[[1]];
sc12[[i,j]]=scatter[[2]];
sc21[[i,j]]=scatter[[3]];
sc22[[i,j]]=scatter[[4]];
\[Kappa]t1[[i,j]]=ff1[[i,j]]+sc11[[i,j]]+sc21[[i,j]];
\[Kappa]t2[[i,j]]=ff2[[i,j]]+sc12[[i,j]]+sc22[[i,j]];
]
]


(* ::Input::Initialization:: *)
For[k=1,k<=10,k++,
For[m=1,m<=Length[\[Tau]],m++,
u1[[m]]=2*\[Pi]/c*(Integrate[Interpolation[Transpose@{\[Theta]x,Sin[\[Theta]x]*(I1u[[m]]+I1d[[m]])},InterpolationOrder->1][t],{t,Min[\[Theta]x],Max[\[Theta]x]}]);
u2[[m]]=2*\[Pi]/c*Integrate[Interpolation[Transpose@{\[Theta]x,Sin[\[Theta]x]*(I2u[[m]]+I2d[[m]])},InterpolationOrder->1][t],{t,Min[\[Theta]x],Max[\[Theta]x]}]
];
(*S1[[m]]=ff1[[m]]/\[Kappa]t1[[m]]*Planck[Es,temperature[\[Tau][[m]]]]/2+sc11[[m]]/\[Kappa]t1[[m]]*(c u1[[m]])/(4 \[Pi])+sc12[[m]]/\[Kappa]t1[[m]]*(c u2[[m]])/(4 \[Pi]);
S2[[m]]=ff2[[m]]/\[Kappa]t2[[m]]*Planck[Es,temperature[\[Tau][[m]]]]/2+sc21[[m]]/\[Kappa]t2[[m]]*(c u1[[m]])/(4 \[Pi])+sc22[[m]]/\[Kappa]t2[[m]]*(c u2[[m]])/(4 \[Pi])];*)

For[j=1,j<=Length[\[Theta]x],j++,
\[Theta]=\[Theta]x[[j]];
\[Mu]=Cos[\[Theta]];
(*upward integration and downward integration*)
For[i=2,i<=Length[\[Tau]],i++,

S1[[i,j]]=ff1[[i,j]]/ \[Kappa]t1[[i,j]]*Planck[Es,temperature[\[Tau][[i]]]]/2+sc11[[i,j]]/\[Kappa]t1[[i,j]]*(c u1[[i]])/(4 \[Pi])+sc12[[i,j]]/\[Kappa]t1[[i,j]]*(c u2[[i]])/(4 \[Pi]);

S2[[i,j]]=ff2[[i,j]]/ \[Kappa]t2[[i,j]]*Planck[Es,temperature[\[Tau][[i]]]]/2+sc21[[i,j]]/\[Kappa]t2[[i,j]]*(c u1[[i]])/(4 \[Pi])+sc22[[i,j]]/\[Kappa]t2[[i,j]]*(c u2[[i]])/(4 \[Pi]);

\[Delta]\[Tau][[i]]=\[Tau][[i]]-\[Tau][[i-1]];
I1u[[i,j]]=1/(1-(\[Delta]\[Tau][[i]]* \[Kappa]t1[[i,j]])/(2*\[Mu]*\[Kappa]T)) ((1+(\[Delta]\[Tau][[i]]* \[Kappa]t1[[i-1,j]])/(2*\[Mu]*\[Kappa]T))I1u[[i-1,j]]-\[Delta]\[Tau][[i]]/(2*\[Mu]*\[Kappa]T)*(\[Kappa]t1[[i-1,j]]*S1[[i-1,j]]+\[Kappa]t1[[i,j]]*S1[[i,j]]));
I2u[[i,j]]=1/(1-(\[Delta]\[Tau][[i]]* \[Kappa]t2[[i,j]])/(2*\[Mu]*\[Kappa]T)) ((1+(\[Delta]\[Tau][[i]]* \[Kappa]t2[[i-1,j]])/(2*\[Mu]*\[Kappa]T))I2u[[i-1,j]]-\[Delta]\[Tau][[i]]/(2*\[Mu]*\[Kappa]T)*(\[Kappa]t2[[i-1,j]]*S2[[i-1,j]]+\[Kappa]t2[[i,j]]*S2[[i,j]]));
ax=I1u[[i,j]];
bx=I2u[[i,j]];

i1=Length[\[Tau]]+1-i;
S1[[i1,j]]=ff1[[i1,j]]/ \[Kappa]t1[[i1,j]]*Planck[Es,temperature[\[Tau][[i1]]]]/2+sc11[[i1,j]]/\[Kappa]t1[[i1,j]]*(c u1[[i1]])/(4 \[Pi])+sc12[[i1,j]]/\[Kappa]t1[[i1,j]]*(c u2[[i1]])/(4 \[Pi]);

S2[[i1,j]]=ff2[[i1,j]]/ \[Kappa]t2[[i1,j]]*Planck[Es,temperature[\[Tau][[i1]]]]/2+sc21[[i1,j]]/\[Kappa]t2[[i1,j]]*(c u1[[i1]])/(4 \[Pi])+sc22[[i1,j]]/\[Kappa]t2[[i1,j]]*(c u2[[i1]])/(4 \[Pi]);
\[Delta]\[Tau]1[[i1]]=\[Tau][[i1]]-\[Tau][[i1+1]];
I1d[[i1,j]]=1/(1+(\[Delta]\[Tau]1[[i1]]* \[Kappa]t1[[i1,j]])/(2*\[Mu]*\[Kappa]T)) ((1-(\[Delta]\[Tau]1[[i1]]* \[Kappa]t1[[i1+1,j]])/(2*\[Mu]*\[Kappa]T))I1d[[i1+1,j]]+\[Delta]\[Tau]1[[i1]]/(2*\[Mu]*\[Kappa]T)*(\[Kappa]t1[[i1+1,j]]*S1[[i1+1,j]]+\[Kappa]t1[[i1,j]]*S1[[i1,j]]));
I2d[[i1,j]]=1/(1+(\[Delta]\[Tau]1[[i1]]* \[Kappa]t2[[i1,j]])/(2*\[Mu]*\[Kappa]T)) ((1-(\[Delta]\[Tau]1[[i1]]* \[Kappa]t2[[i1+1,j]])/(2*\[Mu]*\[Kappa]T))I2d[[i1+1,j]]+\[Delta]\[Tau]1[[i1]]/(2*\[Mu]*\[Kappa]T)*(\[Kappa]t2[[i1+1,j]]*S2[[i1+1,j]]+\[Kappa]t2[[i1,j]]*S2[[i1,j]]));
ay=I1d[[i1,j]];
by=I2d[[i1,j]];

If[i==position,
I1u[[i,j]]=Pjump[Es,B,\[Theta]]*ax+(1-Pjump[Es,B,\[Theta]])*bx;
I2u[[i,j]]=Pjump[Es,B,\[Theta]]*bx+(1-Pjump[Es,B,\[Theta]])*ax;
I1d[[i1,j]]=Pjump[Es,B,\[Theta]]*ay+(1-Pjump[Es,B,\[Theta]])*by;
I2d[[i1,j]]=Pjump[Es,B,\[Theta]]*by+(1-Pjump[Es,B,\[Theta]])*ay,
I1u[[i,j]]=ax;I2u[[i,j]]=bx;I1d[[i1,j]]=ay;I2d[[i1,j]]=by
];
]
];

];


(* ::Input::Initialization:: *)
data=Transpose@{\[Theta]x,I1u[[-1]],I2u[[-1]]};
SetDirectory[NotebookDirectory[]];
Export["third.dat",data];
