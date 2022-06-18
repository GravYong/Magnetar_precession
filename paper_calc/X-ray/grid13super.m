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
\[Tau]mid=27.1;
a1=0.793;
a2= 0.122;
a3=-0.502;
a4=0.548;
a5= -0.205;
a6=0.0266;
b3=0.00445;
b4=0.0108;
b5=0.00211;
b6=0.0000574;
x=Log[10,\[Tau]]-Log[10,\[Tau]mid];
Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]>\[Tau]mid && \[Tau]<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]<= \[Tau]mid && \[Tau]>= 10^-3}}]];
LogLogPlot[temperature[\[Tau]],{\[Tau],0.001,10^4},Background->White,PlotTheme->"Monochrome"]


(* ::Section:: *)
(*Scatter opacity*)


(* ::Input::Initialization:: *)
Aint[Es_,B_,\[Tau]_]:=Module[{Ec,\[Omega],\[Omega]c,Ece,Eci,Epe,Epi,ue,ui,ve,vi,b,a,q,m,ri,\[Epsilon]p,\[Eta]p,g,\[Epsilon],\[Eta],\[Beta]i,K1i,K2i,Kz1i,Kz2i,ep1i,em1i,ep2i,em2i,e01i,e02i,Ap1,Am1,A01,Ap2,Am2,A02,Ap,Am,A0,\[Rho],x,T,Mns,Rns,gns,\[Tau]mid,\[Gamma]eiL,\[Gamma]eip,y1,y2,\[Gamma]re,\[Gamma]ri,gl,gp,\[Epsilon]pg,\[Epsilon]mg,\[Delta]v},
\[Tau]mid=27.1;
a1=0.793;
a2= 0.122;
a3=-0.502;
a4=0.548;
a5= -0.205;
a6=0.0266;
b3=0.00445;
b4=0.0108;
b5=0.00211;
b6=0.0000574;
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
(*turn on the vacuum*)
\[Delta]v=\[Alpha]F/(45 \[Pi]) b^2;
a=1-2*\[Delta]v;
q=7\[Delta]v;
m=-4\[Delta]v;
(*turn off the vacuum*)
(*a=1;
q=0;
m=0;*)
\[Gamma]re=9.5*10^-6*Es/1;
\[Gamma]ri=5.2*10^-9*Z^2/A*Es/1;
y1=hbar*\[Omega]c/(kb*T);
y2=\[Omega]/\[Omega]c;
gl=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*C1[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
gp=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*2*y2*Exp[2x1]C0[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
\[Gamma]eiL=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2)*(1-Exp[-Ec/(kb*T)])*gl;
\[Gamma]eip=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2) (1-Exp[-Ec/(kb*T)])*gp;
\[Epsilon]pg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re+ue^0.5)(1+I*\[Gamma]ri-ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]mg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re-ue^0.5)(1+I*\[Gamma]ri+ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]p=1/2 (\[Epsilon]pg+\[Epsilon]mg);
g=1/2 (\[Epsilon]pg-\[Epsilon]mg);
\[Eta]p=1-ve/(1+I*(\[Gamma]eip+\[Gamma]re))-vi/(1+I*(\[Gamma]eip+\[Gamma]ri));

\[Epsilon]=\[Epsilon]p+a-1;
\[Eta]=\[Eta]p+a+q-1;

ri=1+m/a Sin[\[Theta]i]^2;
\[Beta]i=-((\[Epsilon]^2-g^2-\[Epsilon] *\[Eta](1+m/a))/(2g \[Eta]))*Sin[\[Theta]i]^2/Cos[\[Theta]i];
K1i=\[Beta]i(1-(1+ri/\[Beta]i^2)^(1/2));
K2i=\[Beta]i(1+(1+ri/\[Beta]i^2)^(1/2));
Kz1i=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]i]Cos[\[Theta]i]K1i + g * Sin[\[Theta]i])/(\[Epsilon]p*Sin[\[Theta]i]^2+(\[Eta]p+q)Cos[\[Theta]i]^2+a-1);

Kz2i=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]i]Cos[\[Theta]i]K2i + g * Sin[\[Theta]i])/(\[Epsilon]p*Sin[\[Theta]i]^2+(\[Eta]p+q)Cos[\[Theta]i]^2+a-1);
ep1i=Abs[1+K1i*Cos[\[Theta]i]+Kz1i*Sin[\[Theta]i]]^2/(2(1+Abs[K1i]^2+Abs[Kz1i]^2));
ep2i=Abs[1+K2i*Cos[\[Theta]i]+Kz2i*Sin[\[Theta]i]]^2/(2(1+Abs[K2i]^2+Abs[Kz2i]^2));
em1i=Abs[1-(K1i*Cos[\[Theta]i]+Kz1i*Sin[\[Theta]i])]^2/(2(1+Abs[K1i]^2+Abs[Kz1i]^2));
em2i=Abs[1-(K2i*Cos[\[Theta]i]+Kz2i*Sin[\[Theta]i])]^2/(2(1+Abs[K2i]^2+Abs[Kz2i]^2));
e01i=Abs[K1i * Sin[\[Theta]i]-Kz1i*Cos[\[Theta]i]]^2/(1+Abs[K1i]^2+Abs[Kz1i]^2);
e02i=Abs[K2i * Sin[\[Theta]i]-Kz2i*Cos[\[Theta]i]]^2/(1+Abs[K2i]^2+Abs[Kz2i]^2);
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
scatteropacity[Es_,B_,\[Tau]_,\[Theta]b_]:=Module[{Ec,\[Omega],\[Omega]c,Ece,Eci,Epe,Epi,ue,ui,ve,vi,b,a,q,m,r,\[Epsilon]p,\[Eta]p,g,\[Epsilon],\[Eta],\[Beta],K1,K2,Kz1,Kz2,ep1,em1,ep2,em2,e01,e02,gl,gp,\[Omega]ci,\[Tau]mid,x,T,\[Rho],Mns,Rns,gns,\[Gamma]eiL,\[Gamma]eip,y1,y2,\[Gamma]re,\[Gamma]ri,\[Epsilon]pg,\[Epsilon]mg,\[Gamma]ep,\[Gamma]em,\[Gamma]e0,\[Gamma]ip,\[Gamma]im,\[Gamma]i0,\[Kappa]s1e1,\[Kappa]s1i1,\[Kappa]s2e1,\[Kappa]s2i1,\[Kappa]s11,\[Kappa]s12,\[Kappa]s1e2,\[Kappa]s1i2,\[Kappa]s2e2,\[Kappa]s2i2,\[Kappa]s21,\[Kappa]s22,\[Delta]v},
\[Tau]mid=27.1;
a1=0.793;
a2= 0.122;
a3=-0.502;
a4=0.548;
a5= -0.205;
a6=0.0266;
b3=0.00445;
b4=0.0108;
b5=0.00211;
b6=0.0000574;
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
\[Gamma]re=9.5*10^-6*Es/1;
\[Gamma]ri=5.2*10^-9*Z^2/A*Es/1;
y1=hbar*\[Omega]c/(kb*T);
y2=\[Omega]/\[Omega]c;
gl=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*C1[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
gp=NIntegrate[Exp[-y1*y2*Sinh[x1]^2]*2*y2*Exp[2x1]C0[y2*Exp[2 x1]],{x1,-Infinity,Infinity}];
\[Gamma]eiL=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2)*(1-Exp[-Ec/(kb*T)])*gl;
\[Gamma]eip=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2) (1-Exp[-Ec/(kb*T)])*gp;
\[Epsilon]pg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re+ue^0.5)(1+I*\[Gamma]ri-ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]mg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re-ue^0.5)(1+I*\[Gamma]ri+ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]p=1/2 (\[Epsilon]pg+\[Epsilon]mg);
g=1/2 (\[Epsilon]pg-\[Epsilon]mg);
\[Eta]p=1-ve/(1+I*(\[Gamma]eip+\[Gamma]re))-vi/(1+I*(\[Gamma]eip+\[Gamma]ri));

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
ep1=Abs[1+K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b]]^2/(2(1+Abs[K1]^2+Abs[Kz1]^2));
ep2=Abs[1+K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b]]^2/(2(1+Abs[K2]^2+Abs[Kz2]^2));
em1=Abs[1-(K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b])]^2/(2(1+Abs[K1]^2+Abs[Kz1]^2));
em2=Abs[1-(K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b])]^2/(2(1+Abs[K2]^2+Abs[Kz2]^2));
e01=Abs[K1 * Sin[\[Theta]b]-Kz1*Cos[\[Theta]b]]^2/(1+Abs[K1]^2+Abs[Kz1]^2);
e02=Abs[K2 * Sin[\[Theta]b]-Kz2*Cos[\[Theta]b]]^2/(1+Abs[K2]^2+Abs[Kz2]^2);
\[Gamma]ep=\[Gamma]eiL+\[Gamma]re;
\[Gamma]em=\[Gamma]eiL+\[Gamma]re;
\[Gamma]e0=\[Gamma]eip+\[Gamma]re;
\[Gamma]ip=me/mp (\[Gamma]eiL+\[Gamma]ri);
\[Gamma]im=me/mp (\[Gamma]eiL+\[Gamma]ri);
\[Gamma]i0=me/mp (\[Gamma]eip+\[Gamma]ri);
\[Kappa]s1e1=Ye*\[Kappa]T*(((1+ue^0.5)^2+\[Gamma]ep^2)^-1 ep1*Ap1+((1-ue^0.5)^2+\[Gamma]em^2)^-1 em1*Am1+((1+0*ue^0.5)^2+\[Gamma]e0^2)^-1 e01*A01);
\[Kappa]s1i1=((Z^2 me)/(A mp))^2*\[Kappa]T/A*(((1-ui^0.5)^2+\[Gamma]ip^2)^-1 ep1*Ap1+((1+ui^0.5)^2+\[Gamma]im^2)^-1 em1*Am1+((1+0*ui^0.5)^2+\[Gamma]i0^2)^-1 e01*A01);
\[Kappa]s11=\[Kappa]s1e1+\[Kappa]s1i1;

\[Kappa]s1e2=Ye*\[Kappa]T*(((1+ue^0.5)^2+\[Gamma]ep^2)^-1 ep1*Ap2+((1-ue^0.5)^2+\[Gamma]em^2)^-1 em1*Am2+((1+0*ue^0.5)^2+\[Gamma]e0^2)^-1 e01*A02);
\[Kappa]s1i2=((Z^2 me)/(A mp))^2*\[Kappa]T/A*(((1-ui^0.5)^2+\[Gamma]ip^2)^-1 ep1*Ap2+((1+ui^0.5)^2+\[Gamma]im^2)^-1 em1*Am2+((1+0*ui^0.5)^2+\[Gamma]i0^2)^-1 e01*A02);
\[Kappa]s12=\[Kappa]s1e2+\[Kappa]s1i2;

\[Kappa]s2e1=Ye*\[Kappa]T*(((1+ue^0.5)^2+\[Gamma]ep^2)^-1 ep2*Ap1+((1-ue^0.5)^2+\[Gamma]em^2)^-1 em2*Am1+((1+0*ue^0.5)^2+\[Gamma]e0^2)^-1 e02*A01);
\[Kappa]s2i1=((Z^2 me)/(A mp))^2*\[Kappa]T/A*(((1-ui^0.5)^2+\[Gamma]ip^2)^-1 ep2*Ap1+((1+ui^0.5)^2+\[Gamma]im^2)^-1 em2*Am1+((1+0*ui^0.5)^2+\[Gamma]i0^2)^-1 e02*A01);
\[Kappa]s21=\[Kappa]s2e1+\[Kappa]s2i1;

\[Kappa]s2e2=Ye*\[Kappa]T*(((1+ue^0.5)^2+\[Gamma]ep^2)^-1 ep2*Ap2+((1-ue^0.5)^2+\[Gamma]em^2)^-1 em2*Am2+((1+0*ue^0.5)^2+\[Gamma]e0^2)^-1 e02*A02);
\[Kappa]s2i2=((Z^2 me)/(A mp))^2*\[Kappa]T/A*(((1-ui^0.5)^2+\[Gamma]ip^2)^-1 ep2*Ap2+((1+ui^0.5)^2+\[Gamma]im^2)^-1 em2*Am2+((1+0*ui^0.5)^2+\[Gamma]i0^2)^-1 e02*A02);
\[Kappa]s22=\[Kappa]s2e2+\[Kappa]s2i2;
{\[Kappa]s11,\[Kappa]s12,\[Kappa]s21,\[Kappa]s22}
]


(* ::Section:: *)
(*free free opacity*)


(* ::Input::Initialization:: *)
freeopacity[Es_,B_,\[Tau]_,\[Theta]b_]:=Module[{\[Epsilon]p,\[Eta]p,\[Epsilon],\[Eta],g,Ec,\[Omega],\[Omega]c,Ece,Epi,Eci,Epe,ue,ui,ve,vi,\[Gamma]re,\[Gamma]ri,r1,r2,\[Gamma]eiL,\[Gamma]eip,\[Epsilon]pg,\[Epsilon]mg,\[Alpha]F,BQ,b,a,q,m,r,\[Beta],K1,K2,Kz1,Kz2,ep1,ep2,em1,em2,e01,e02,\[Gamma]p,\[Gamma]m,\[CapitalLambda]p,\[CapitalLambda]m,\[Kappa]p,\[Kappa]m,\[Kappa]o,\[Kappa]ff1,\[Kappa]ff2,\[Tau]mid,x,Mns,Rns,gns,gp,gl,\[Rho],T,\[Delta]v},
\[Tau]mid=27.1;
a1=0.793;
a2= 0.122;
a3=-0.502;
a4=0.548;
a5= -0.205;
a6=0.0266;
b3=0.00445;
b4=0.0108;
b5=0.00211;
b6=0.0000574;
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
\[Gamma]re=9.5*10^-6*Es/1;
\[Gamma]ri=5.2*10^-9*Z^2/A*Es/1;
r1=hbar*\[Omega]c/(kb*T);
r2=\[Omega]/\[Omega]c;
gl=NIntegrate[Exp[-r1*r2*Sinh[x1]^2]*C1[r2*Exp[2 x1]],{x1,-Infinity,Infinity}];
gp=NIntegrate[Exp[-r1*r2*Sinh[x1]^2]*2*r2*Exp[2x1]C0[r2*Exp[2 x1]],{x1,-Infinity,Infinity}];
\[Gamma]eiL=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2)*(1-Exp[-Ec/(kb*T)])*gl;
\[Gamma]eip=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2) (1-Exp[-Ec/(kb*T)])*gp;

\[Epsilon]pg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re+ue^0.5)(1+I*\[Gamma]ri-ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]mg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re-ue^0.5)(1+I*\[Gamma]ri+ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]p=1/2 (\[Epsilon]pg+\[Epsilon]mg);
g=1/2 (\[Epsilon]pg-\[Epsilon]mg);
\[Eta]p=1-ve/(1+I*(\[Gamma]eip+\[Gamma]re))-vi/(1+I*(\[Gamma]eip+\[Gamma]ri));

\[Alpha]F=e^2/hbar/c;
BQ=me^2*c^3/e/hbar;
\[Epsilon]=\[Epsilon]p+a-1;
\[Eta]=\[Eta]p+a+q-1;
r=1+m/a Sin[\[Theta]b]^2;
\[Beta]=-((\[Epsilon]^2-g^2-\[Epsilon] *\[Eta](1+m/a))/(2g \[Eta]))*Sin[\[Theta]b]^2/Cos[\[Theta]b];
K1=\[Beta](1-(1+r/\[Beta]^2)^(1/2));
K2=\[Beta](1+(1+r/\[Beta]^2)^(1/2));
Kz1=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K1 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1);

Kz2=-((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K2 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1);
ep1=Abs[1+K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b]]^2/(2(1+Abs[K1]^2+Abs[Kz1]^2));
ep2=Abs[1+K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b]]^2/(2(1+Abs[K2]^2+Abs[Kz2]^2));
em1=Abs[1-(K1*Cos[\[Theta]b]+Kz1*Sin[\[Theta]b])]^2/(2(1+Abs[K1]^2+Abs[Kz1]^2));
em2=Abs[1-(K2*Cos[\[Theta]b]+Kz2*Sin[\[Theta]b])]^2/(2(1+Abs[K2]^2+Abs[Kz2]^2));
e01=Abs[K1 * Sin[\[Theta]b]-Kz1*Cos[\[Theta]b]]^2/(1+Abs[K1]^2+Abs[Kz1]^2);
e02=Abs[K2 * Sin[\[Theta]b]-Kz2*Cos[\[Theta]b]]^2/(1+Abs[K2]^2+Abs[Kz2]^2);

\[Gamma]p=\[Gamma]eiL+(1+ue^0.5)\[Gamma]ri+(1-ui^0.5)\[Gamma]re;
\[Gamma]m=\[Gamma]eiL+(1-ue^0.5)\[Gamma]ri+(1+ui^0.5)\[Gamma]re;
\[CapitalLambda]p=((1+ue^0.5)^2 (1-ui^0.5)^2+\[Gamma]p^2)^-1;
\[CapitalLambda]m=((1-ue^0.5)^2 (1+ui^0.5)^2+\[Gamma]m^2)^-1;
\[Kappa]p=\[Omega]/( c \[Rho]) ve \[CapitalLambda]p \[Gamma]eiL;
\[Kappa]m=\[Omega]/( c \[Rho]) ve \[CapitalLambda]m \[Gamma]eiL;
\[Kappa]o=\[Omega]/( c \[Rho]) ve*\[Gamma]eip;
\[Kappa]ff1=\[Kappa]p*ep1 + \[Kappa]m*em1 + \[Kappa]o*e01;
\[Kappa]ff2=\[Kappa]p*ep2 + \[Kappa]m*em2+ \[Kappa]o*e02;
{\[Kappa]ff1,\[Kappa]ff2}]


(* ::Input::Initialization:: *)
freeopacity1[Es_,B_,\[Tau]_,\[Theta]b_]:=Module[{Ec,\[Omega],\[Omega]c,Ece,Eci,Epe,Epi,ue,ui,ve,vi,b,a,q,m,r,\[Epsilon]p,\[Eta]p,g,\[Epsilon],\[Eta],\[Beta],K1,K2,Kz1,Kz2,ep1,em1,ep2,em2,e01,e02,ne,ni,\[Alpha]0,gl,gp,y1,y2,\[Omega]ci,\[Nu]re,\[Nu]ep,\[Nu]em,\[Nu]e0,\[Nu]ip,\[Nu]im,\[Nu]i0,ffe1,ffi1,ffe2,ffi2,ff1,ff2,T,\[Rho],Mns,Rns,gns,\[Tau]mid,x,\[Delta]v},
\[Tau]mid=27.1;
a1=0.793;
a2= 0.122;
a3=-0.502;
a4=0.548;
a5= -0.205;
a6=0.0266;
b3=0.00445;
b4=0.0108;
b5=0.00211;
b6=0.0000574;
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
\[Delta]v=\[Alpha]F/(45 \[Pi]) b^2;
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
\[Delta]v=\[Alpha]F/(45 \[Pi]) b^2;
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
Es=1;
B=10^13;
(*\[Tau]v=\[Tau]reson[Es,B];
\[Tau]x1=logspace[-3,Log[10,\[Tau]v],(IntegerPart[Log[10,\[Tau]v]]+3)*20];
\[Tau]x2=logspace[Log[10,\[Tau]v],Log[10,2]+4,IntegerPart[(Log[10,2]+4-Log[10,\[Tau]v])]*20];
\[Tau]=Reverse[DeleteDuplicates[Flatten[{\[Tau]x1,\[Tau]x2}]]];
position=IntegerPart[(Log[10,2]+4-Log[10,\[Tau]v])]*20;*)
\[Tau]=Reverse[logspace[-3,4,2000]];
\[Theta]x=N[Range[0.0001,89.9,(89.9-0.0001)/30]]/180*\[Pi];
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
(*For[j=1,j\[LessEqual]Length[\[Theta]x],j++,
For[
i=1,i\[LessEqual]Length[\[Tau]],i++,S1[[i,j]]=1/2*Planck[Es,temperature[\[Tau][[i]]]];S2[[i,j]]=1/2*Planck[Es,temperature[\[Tau][[i]]]]];
]*)


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
scatter=scatteropacity[Es,B,\[Tau][[i]],\[Theta]x[[j]]];
sc11[[i,j]]=scatter[[1]];
sc12[[i,j]]=scatter[[2]];
sc21[[i,j]]=scatter[[3]];
sc22[[i,j]]=scatter[[4]];
\[Kappa]t1[[i,j]]=ff1[[i,j]]+sc11[[i,j]]+sc12[[i,j]];
\[Kappa]t2[[i,j]]=ff2[[i,j]]+sc21[[i,j]]+sc22[[i,j]];
]
]


(* ::Input::Initialization:: *)
Export["ff1.dat",ff1];
Export["ff2.dat",ff2];
Export["sc11.dat",sc11];
Export["sc12.dat",sc12];
Export["sc21.dat",sc21];
Export["sc22.dat",sc22];
