(* ::Package:: *)

(* ::Input:: *)
(*(*Values*)*)
(*(*ehats=0.6;*)
(*ehato=0.6;*)
(*ws=1;*)
(*ms=1;*)
(*ps=5000;*)
(*qs=0.002;*)
(*wo=1;*)
(*mo=1;*)
(*po=5000;*)
(*qo=0.002;*)
(*a=-0.03;*)
(*b=0.1102;*)
(*eco=0;*)
(*crit=0.1;*)
(*k=10;*)
(*r=0.2;*)
(*gammao=0.5;*)
(*gammas=0.6;*)*)
(**)
(*valueList = {gammao->0.3,gammas->0.5,ehats->0.6,ehato-> 0.6,ws-> 1,ms-> 1,ps-> 5000,qs-> 0.002,wo-> 1,mo-> 1,po-> 5000,qo-> 0.002,a-> -0.03,b-> 0.1102,eco-> 0,crit-> 0.1,k-> 10,r-> 0.2};*)
(*trueVanilla = {r->0.2,k->10,crit->0.1,ns->100,no->100,ps->5000,po->5000,qs->0.007,qo->0.007,ws->1,wo->1,ms->1,mo->1,ehats->0.6,ehato->0.6,gammas->0.1,gammao->0.1,betas->0.2,betas->0.2,a->-0.3,b->0.1102};*)
(**)
(*(*Expressions*)*)
(*ns:=(Cs+ds);*)
(*no :=(Co+do);*)
(**)
(*(*Functions of (x,t,Cs,Co)*)*)
(*Rs[x_]:=a+b*x;*)
(*hcs[x_]:=qs*x*ecs[x];*)
(*hco[x_]:=qs*x*eco;*)
(*hds[x_]:=qs*x*eds[x];*)
(*hdo[x_]:=qs*x*edo[x];*)
(*eds[x_]:= Piecewise[{{ehats,x >= (ws+ms)/(ps*qs)},{0,x<(ws+ms)/(ps*qs)}}];*)
(*ecs[x_]:=Max[0,{Min[ehats,Rs[x]]}];*)
(*edo[x_]:= Piecewise[{{ehato,x >= (wo+mo)/(po*qo)},{0,x <(wo+mo)/(po*qo)}}];*)
(*pics[x_]:=ps*hcs[x]-ws*ecs[x]+ms*(ehats-ecs[x]);*)
(*pids[x_]:=ps*hds[x]-ws*eds[x]+ms*(ehats-eds[x]);*)
(*pico[x_]:=po*hco[x]-wo*eco+mo*(ehato-eco);*)
(*pido[x_]:=po*hdo[x]-wo*edo[x]+mo*(ehato-edo[x]);*)
(*Eff[x_,Cs_,Co_]:=qs*(eds[x]*(ns-Cs)+qs*ecs[x]*Cs)+qo*(edo[x]*(no-Co)+qo*eco*Cs);*)
(**)
(*(*Reduced Functions for X's Fixed Points*)*)
(*zeta[x_,Cs_,Co_] :=Eff[x,Cs,Co]-qs*Rs*Cs;*)
(*kappa := -r/k;*)
(*omega[Cs_]:=r*crit+r-b*qs*Cs;*)
(*chi[x_,Cs_,Co_]:=-r*crit-a*qs*Cs-zeta[x,Cs,Co];*)
(**)
(*(*Fixed Points*)*)
(*xfp2[x_,Cs_,Co_]:=(omega[Cs]+(omega[Cs]^2-4*kappa*chi[x,Cs,Co])^(1/2))(2*kappa)*)
(*csfp2[x_]:=ns*((-betas/gammas)*(1-(pics[x]/pids[x]))+1);*)
(*cofp2[x_]:=no*((-betao/gammao)*(1-(pico[x]/pido[x]))+1);*)
(**)
(*(*Partial Dervatives for use in Jacobian calculated by hand.*)*)
(*(*df1dx is in a pre-reduced form, it should technically have Deff as well, but that can be computed pre-solution.*)*)
(*df1dx[x_]:=r*x-(3*r*x^3)/k-r*crit+(2*r*crit*x)/k-Eff[x,Cs,Co]*)
(*df1dCs[x_]:=qs*x*(eds[x]-ecs[x]);*)
(*df1dCo[x_]:=qo*x*(edo[x]-eco);*)
(*df2dx[x_,Cs_]:=betas*Cs*(Dpics[x]*pids[x]-pics[x]*Dpids[x])/(pids[x]^2);*)
(*df2dCs[x_,Cs_]:=gammas-(2*gammas*Cs)/ns*betas*(1-(pics[x]/pids[x]));*)
(*df3dx[x_,Co_]:=betao*Co*(Dpico[x]*pido[x]-pico[x]*Dpido[x])/(pido[x]^2);*)
(*df3dCo[x_,Co_]:=gammao-(2*gammao*Co)/no*betao*(1-(pico[x]/pido[x]));*)
(**)
(*(*Fixed Point Substitution Sets*)*)
(*fpSet1 = {x->0,Cs->0,Co->0};*)
(*fpSet2={x->0,Cs->csfp2[x]/.{x->0},Co->0};*)
(*fpSet3={x->0 ,Cs->0, Co->cofp2[x]/.{x->0}};*)
(*fpSet4={x->0,Cs->csfp2[x]/.{x->0},Co->cofp2[x]/.{x->0}};*)
(*fpSet5={x->xTWO,Cs->0,Co->0};*)
(*fpSet6={x->xTWO,Cs->CSTWO,Co->0};*)
(*fpSet7={x->xTWO,Cs->0,Co->COTWO};*)
(**)
(*(*Embed partial DE's into matrix form.*)*)
(*symbolicMatrix = {{df1dx[x],df1dCs[x],df1dCo[x]},{df2dx[x,Cs],df2dCs[x,Cs],0},{df3dx[x,Co],0,df3dCo[x,Co]}};*)
(*displayMatrix = MatrixForm[{{df1dx,df1dCs,df1dCo},{df2dx,df2dCs,0},{df3dx,0,df3dCo}}]*)
(*(*numericMatrix = symbolicMatrix/.Symbol["fpSet"<>ToString[7]];*)*)
(*MatrixForm[FullSimplify[Eigenvalues[numericMatrix]]];*)
(**)
(*(*Loops and solves all fixed point eigenvalues. Warning: The last set (ie. fpset7) is VERY computationally expensive, could take a few hours.*)*)
(*Symbol["fpSet"<>ToString[1]]*)
(*For[i=1,i<=6,i++,Print["FPSet ",i, " ::->::",MatrixForm[Symbol["fpSet"<>ToString[i]]],"Eival ::->::",MatrixForm[FullSimplify[Eigenvalues[symbolicMatrix/.Symbol["fpSet"<>ToString[i]]]]],"\n\n\n"]]*)
(**)
(**)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
