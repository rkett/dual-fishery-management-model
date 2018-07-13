(* ::Package:: *)

(* ::Input:: *)
(*(* This is a special workbook dedicated to solving just x2,x3,x4,x5 fixed points with cs1,cs2,co1 and co2. I've done it like this as my original workbook for fixed point solutions of a jacobian determinant was a bit too messy, and,from by-hand analysis I can make pre-reductions to certain equations before any substitution is done.*)
(**)*)
(**)
(*trueVanilla = {r->0.2,k->10,crit->0.1,ns->100,no->100,ps->5000,po->5000,qs->0.007,qo->0.007,ws->1,wo->1,ms->1,mo->1,ehats->0.6,ehato->0.6,gammas->0.1,gammao->0.1,betas->0.2,betao->0.2,a->-0.3,b->0.1102};*)
(**)
(*(*Special Pre-Reduced form of this for x4,x5; as: Deff \[Equal] 0 for all R_S > ehats (ie. x4,x5 have to fufill this condition*)*)
(*(*df1dx[x_]:=2*r*x-(3*r*x^2)/k-r*crit+(2*r*crit*x)/k-Eff[x,Cs,Co]; *)*)
(**)
(*(*Special Pre-Reduced form of this for x2,x3; as: Deff \[Equal] qsbCs for all 0 < R_S < ehats (ie. x2,x3 have to fufill this condition*)*)
(*df1dx[x_]:=2*r*x-(3*r*x^2)/k-r*crit+(2*r*crit*x)/k-Eff[x,Cs,Co]-qs*b*Cs;*)
(**)
(*(*Normal When Compared to General DetJ file.*)*)
(*df1dCs[x_]:=qs*x*(eds[x]-ecs[x]);*)
(*df1dCo[x_]:=qo*x*(edo[x]-eco);*)
(*df2dx[x_,Cs_]:=betas*Cs*(Dpics[x]*pids[x]-pics[x]*Dpids[x])/(pids[x]^2);*)
(*df2dCs[x_,Cs_]:=gammas-(2*gammas*Cs)/ns*betas*(1-(pics[x]/pids[x]));*)
(*df3dx[x_,Co_]:=betao*Co*(Dpico[x]*pido[x]-pico[x]*Dpido[x])/(pido[x]^2);*)
(*df3dCo[x_,Co_]:=gammao-(2*gammao*Co)/no*betao*(1-(pico[x]/pido[x]));*)
(**)
(*(*zeta[x_,Cs_,Co_] :=Eff[x,Cs,Co]-qs*ehats*Cs;*)
(*kappa := -r/k;*)
(*tau :=(r*crit)/k+r;*)
(*rho[Cs_]:=-r*crit-ehats*qs*Cs-zeta[x,Cs,Co];*)*)
(**)
(*xfp4[Cs_,Co_]:=-(tau+(tau^2-4*kappa*rho[Cs])^(1/2))/(2*kappa);*)
(*xfp5[x_,Cs_,Co_]:=-(tau-(tau^2-4*kappa*rho[Cs])^(1/2))/(2*kappa);*)
(**)
(*(*CS2=ns*(-(betas/gammas)*(1-pics[x]/pids[x])+1);*)*)
(**)
(*fpSet5 = {x->x/.{Cs->CS2,Co->0},Cs->CS2,Co->0};*)
(*fpSet6 = {x->x/.{Cs->0,Co->CO2},Cs->0,Co->CO2};*)
(*fpSet7= {x->x/.{Cs->CS2,Co->CO2},Cs->CS2,Co->CO2};*)
(**)
(*symbolicMatrix = {{df1dx[x],df1dCs[x],df1dCo[x]},{df2dx[x,Cs],df2dCs[x,Cs],0},{df3dx[x,Co],0,df3dCo[x,Co]}};*)
(*displayMatrix = MatrixForm[{{df1dx,df1dCs,df1dCo},{df2dx,df2dCs,0},{df3dx,0,df3dCo}}];*)
(*numericMatrix = symbolicMatrix/.Symbol["fpSet"<>ToString[6]];*)
(**)
(*EIGENVALUES =MatrixForm[FullSimplify[Eigenvalues[numericMatrix]]]*)
(*(*EIGENVECTORS =MatrixForm[FullSimplify[Eigenvectors[numericMatrix]]]*)*)
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
