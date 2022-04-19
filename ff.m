(* ::Package:: *)

(* ::Input::Initialization:: *)
hbar=1.0545718176461565*10^(-27);
e=4.803204712570263*10^(-10);
me=9.1093837015*10^(-28);
c=29979245800.0;
kb=1.380649*10^(-16);
kev2erg=1.6021766339999998*10^-9;
mp=1.67262192369*10^-24;
G=6.674299999999999*10^-8;
Ms=1.988409870698051*10^33;
\[Sigma]T=6.65*10^-25;
\[Kappa]T=\[Sigma]T/mp;
logspace[a_,b_,n_]:=10.0^Range[a,b,(b-a)/(n-1)];
linearspace[a_,b_,n_]:=Range[a,b,(b-a)/(n-1)];
Planck[Es_,T_]:=2*hbar*2*\[Pi]*(Es*kev2erg/hbar/2/\[Pi])^3/c^2*1/(Exp[(Es*kev2erg)/(kb*T)]-1);

C0[a_]:=N[1/a-Exp[a]*ExpIntegralE[1,a]]
C1[a_]:=N[(1+a)*Exp[a]*ExpIntegralE[1,a]-1]
freeopacity[B_,Es_,\[Tau]_,\[Theta]b_]:=Module[{\[Epsilon]p,\[Eta]p,\[Epsilon],\[Eta],g,Z,A,Ye,M,Ec,\[Omega],\[Omega]c,Ece,Epi,Eci,Epe,ue,ui,ve,vi,\[Gamma]re,\[Gamma]ri,r1,r2,\[Gamma]eiL,\[Gamma]eip,\[Epsilon]pg,\[Epsilon]mg,\[Alpha]F,BQ,b,a,q,m,r,\[Beta],K1,K2,Kz1,Kz2,ep1,ep2,em1,em2,e01,e02,\[Gamma]p,\[Gamma]m,\[CapitalLambda]p,\[CapitalLambda]m,\[Kappa]p,\[Kappa]m,\[Kappa]o,\[Kappa]ff1,\[Kappa]ff2,\[Tau]mid,x,Mns,Rns,gns,\[Rho],T,gp,gl},

Z=1;A=1;Ye=Z/A;M=(Z me)/(A mp);
\[Tau]mid=63.2;
a1=0.761;
a2= 0.00198;
a3=0.267;
a4=\[Minus]0.356;
a5= 0.179;
a6=\[Minus]0.0282;
b3=\[Minus]0.118;
b4=\[Minus]0.0336;
b5=\[Minus]0.00428;
b6=\[Minus]0.00022;
x=Log[10,\[Tau]]-Log[10,\[Tau]mid];
T=Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]>\[Tau]mid && \[Tau]<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]<= \[Tau]mid && \[Tau]>= 10^-3}}];
Mns=1.4*Ms;
Rns=10^6;
gns=(G*Mns)/Rns^2 (1-(2G Mns)/(Rns c^2))^(-1/2);
\[Rho]=mp*gns/(2*\[Kappa]T*kb*T)*\[Tau];

Ec=Es*kev2erg;
\[Omega]=Ec/hbar;
\[Omega]c=(e B)/(me c);
Ece=hbar*\[Omega]c/kev2erg;Eci=Ece*M;Epe=hbar*(4*\[Pi]*\[Rho]/mp * Ye * e^2/me)^0.5/kev2erg;
Epi=M^0.5*Epe;
ue=(Ece/Es)^2;ui=(Eci/Es)^2;ve=(Epe/Es)^2;vi=(Epi/Es)^2;
\[Gamma]re=9.5*10^-6*Es/1;
\[Gamma]ri=5.2*10^-9*Z^2/A*Es/1;
r1=hbar*\[Omega]c/(kb*T);
r2=\[Omega]/\[Omega]c;
gl=NIntegrate[Exp[-r1*r2*Sinh[x]^2]*C1[r2*Exp[2 x]],{x,-Infinity,Infinity}];
gp=NIntegrate[Exp[-r1*r2*Sinh[x]^2]*2*r2*Exp[2x]C0[r2*Exp[2 x]],{x,-Infinity,Infinity}];
\[Gamma]eiL=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2)*(1-Exp[-Ec/(kb*T)])*gl;
\[Gamma]eip=9.24*10^-5*(Z^2*(\[Rho]/1))/(A (T/10^6)^0.5 Es^2) (1-Exp[-Ec/(kb*T)])*gp;

\[Epsilon]pg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re+ue^0.5)(1+I*\[Gamma]ri-ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]mg=1-(ve*(1+I * \[Gamma]ri)+vi*(1+I*\[Gamma]re))/((1+I*\[Gamma]re-ue^0.5)(1+I*\[Gamma]ri+ui^0.5)+I*\[Gamma]eiL);
\[Epsilon]p=1/2 (\[Epsilon]pg+\[Epsilon]mg);
g=1/2 (\[Epsilon]pg-\[Epsilon]mg);
\[Eta]p=1-ve/(1+I*(\[Gamma]eip+\[Gamma]re))-vi/(1+I*(\[Gamma]eip+\[Gamma]ri));

\[Alpha]F=e^2/hbar/c;
BQ=me^2*c^3/e/hbar;
b=B/BQ;
a=1+\[Alpha]F/(2\[Pi]) (1.195-2/3 Log[b]-1/b (0.8553+Log[b])-1/(2b^2));
q=-(\[Alpha]F/(2\[Pi]))(-(2/3)b+1.272-1/b (0.3070+Log[b])-0.7003 1/b^2);
m=-(\[Alpha]F/(2\[Pi]))(2/3+1/b (0.1447-Log[b])-1/b^2);
\[Epsilon]=\[Epsilon]p+a-1;
\[Eta]=\[Eta]p+a+q-1;
r=1+m/a Sin[\[Theta]b]^2;
\[Beta]=-((\[Epsilon]^2-g^2-\[Epsilon] *\[Eta](1+m/a))/(2g \[Eta]))*Sin[\[Theta]b]^2/Cos[\[Theta]b];
K1=\[Beta](1-(1+r/\[Beta]^2)^(1/2));
K2=\[Beta](1+(1+r/\[Beta]^2)^(1/2));
Kz1=-(((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K1 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1));

Kz2=-(((\[Epsilon]p-\[Eta]p-q)*Sin[\[Theta]b]Cos[\[Theta]b]K2 + g * Sin[\[Theta]b])/(\[Epsilon]p*Sin[\[Theta]b]^2+(\[Eta]p+q)Cos[\[Theta]b]^2+a-1));
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
\[Kappa]p=Ec/(hbar c \[Rho]) ve \[CapitalLambda]p \[Gamma]eiL;
\[Kappa]m=Ec/(hbar c \[Rho]) ve \[CapitalLambda]m \[Gamma]eiL;
\[Kappa]o=Ec/(hbar c \[Rho]) ve \[Gamma]eip;

\[Kappa]ff1=\[Kappa]p*ep1 + \[Kappa]m*em1 + \[Kappa]o*e01;
\[Kappa]ff2=\[Kappa]p*ep2 + \[Kappa]m*em2+ \[Kappa]o*e02;
{\[Kappa]ff1,\[Kappa]ff2}]

temperature[\[Tau]1_]:=Module[{\[Tau]mid,x},
\[Tau]mid=63.2;
a1=0.761;
a2= 0.00198;
a3=0.267;
a4=\[Minus]0.356;
a5= 0.179;
a6=\[Minus]0.0282;
b3=\[Minus]0.118;
b4=\[Minus]0.0336;
b5=\[Minus]0.00428;
b6=\[Minus]0.00022;
x=Log[10,\[Tau]1]-Log[10,\[Tau]mid];
Piecewise[{{N[10^6*10^(a1+a2*x+a3*x^2+a4*x^3+a5*x^4+a6*x^5)],\[Tau]1>\[Tau]mid && \[Tau]1<=2*10^4},{N[10^6*10^(a1+a2*x+b3*x^2+b4*x^3+b5*x^4+b6*x^5)],\[Tau]1<= \[Tau]mid && \[Tau]1>= 10^-3}}]];

\[Tau]reson[Es_,B_]:=Module[{q,m,\[Alpha]F,\[Omega]c,Ebi,M,R,Ead,Z,A,Ye,BQ,b,\[Delta]v,fB,g,\[Rho]v,\[Tau]v},
Z=1;
A=1;
Ye=Z/A;
\[Alpha]F=e^2/hbar/c;
BQ=me^2*c^3/e/hbar;
b=B/BQ;
q=-(\[Alpha]F/(2\[Pi]))(-(2/3)b+1.272-1/b (0.3070+Log[b])-0.7003 1/b^2);
m=-(\[Alpha]F/(2\[Pi]))(2/3+1/b (0.1447-Log[b])-1/b^2);
\[Delta]v=\[Alpha]F/(45*\[Pi]) b^2;
fB=((3*\[Delta]v)/(q+m))^0.5;
M=1.4*Ms;
R=10^6;
g=(G*M)/R^2 (1-(2G M)/(R c^2))^(-1/2);
\[Omega]c=(e B)/(me c);
Ebi=hbar*\[Omega]c*(me/mp)/kev2erg;
\[Rho]v=0.96* Ye^-1 (Es/1)^2 (B/10^14)^2 fB^-2;
\[Tau]v=\[Tau]1/.FindRoot[temperature[\[Tau]1]-\[Tau]1*(mp g)/(\[Kappa]T 2 \[Rho]v kb)==0 ,{\[Tau]1,((mp g)/(\[Kappa]T 2 \[Rho]v kb*5*10^6))^-1}];
\[Tau]v]

Es=2;
B=5*10^14;
\[Tau]v=\[Tau]reson[Es,B];
logspace[a_,b_,n_]:=10.0^Range[a,b,(b-a)/(n-1)];
linearspace[a_,b_,n_]:=Range[a,b,(b-a)/(n-1)];
\[Tau]x1=logspace[-3,Log[10,\[Tau]v],(IntegerPart[Log[10,\[Tau]v]]+3)*20];
\[Tau]x2=logspace[Log[10,\[Tau]v],4,(4-IntegerPart[Log[10,\[Tau]v]])*20];
\[Tau]=DeleteDuplicates[Flatten[{\[Tau]x1,\[Tau]x2}]];
\[Theta]x=N[Range[1,89,(89-1)/20]]/180*\[Pi];
ff1=Table[1,{Length[\[Tau]]},{Length[\[Theta]x]}];
ff2=Table[1,{Length[\[Tau]]},{Length[\[Theta]x]}];

For[k=1,k<=Length[\[Theta]x],k++,
For[j=1,j<=Length[\[Tau]],j++,free=freeopacity[B,Es,\[Tau][[j]],\[Theta]x[[k]]];
ff1[[j,k]]=free[[1]];ff2[[j,k]]=free[[2]];
]
]
(*SetDirectory[NotebookDirectory[]];*)
Export["ff1.dat",ff1];
Export["ff2.dat",ff2];
