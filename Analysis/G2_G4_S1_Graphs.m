(* ::Package:: *)

(* ::Input:: *)
(*(*Graph of  (\[NoBreak]0.1*)
(*0.1*)
(*1/10 (-0.2 (1.\[VeryThinSpace]-20.2 x+3 x^2)-7 edo[x]-7 eds[x])*)
(**)
(*\[NoBreak]) for trueVanilla values.*)*)
(*(*E for reference: Eff[x_,Cs_,Co_]:=qs*(eds[x]*(ns-Cs)+qs*ecs[x]*Cs)+qo*(edo[x]*(no-Co)+qo*eco*Cs)*)*)
(**)
(*trueVanilla = {r->0.2,k->10,crit->0.1,ns->100,no->100,ps->5000,po->5000,qs->0.007,qo->0.007,ws->1,wo->1,ms->1,mo->1,ehats->0.6,ehato->0.6,gammas->0.1,gammao->0.1,betas->0.2,betao->0.2,a->-0.3,b->0.1102};*)
(**)
(*(*Graph under 0 < R_S < ehats constraint.*) *)
(*(*Pre analysis*)*)
(*(*Min x \[Equal] 2.722323049001815*)*)
(*-a/b/.trueVanilla;*)
(*(*Max x \[Equal] 8.166969147005444*)*)
(*(ehats-a)/b/.trueVanilla;*)
(*(*Cond for eds \[Equal] 0.05714285714285715*)*)
(*(ws+ms)/(ps*qs)/.trueVanilla;*)
(*(*Cond for edo \[Equal] 0.05714285714285715*)*)
(*(wo+mo)/(po*qo)/.trueVanilla;*)
(*(*Conclusion for pre-analysis, eds = edo = max for this constraint*)*)
(*(*Then, to graph this, only requirement is to scale x(t) for EV3!*)*)
(*Plot[{1/10*(-0.2(1-20.2*x+3*x^2)-7*ehato-7*ehats)}/.trueVanilla,{x,(-a/b), (ehats-a)/b}/.trueVanilla, AxesLabel->{Subscript[SuperStar[x], 2][t],Subscript[\[Lambda], 3]}, PlotLabel->"Graph of \!\(\*SubscriptBox[\(\[Lambda]\), \(3\)]\)(\!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t),0,0) versus Real \!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t)"]*)
(**)
(*(*Graph under R_S > ehats constraint.*)*)
(*Plot[{1/10*(-0.2(1-20.2*x+3*x^2)-7*ehato-7*ehats)}/.trueVanilla,{x,(ehats-a)/b, k}/.trueVanilla, AxesLabel->{Subscript[SuperStar[x], 4][t],Subscript[\[Lambda], 3]}, PlotLabel->"Graph of \!\(\*SubscriptBox[\(\[Lambda]\), \(3\)]\)(\!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(4\)]\)(t),0,0) versus Real \!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(4\)]\)(t)"]*)
(**)
(**)


(* ::Input:: *)
(*Plot[{1/10 (-0.2 (1-20.2 x+3 x^2)-7 ehato-7 ehats)}/. trueVanilla,{x,-(a/b),(ehats-a)/b}/. trueVanilla,PlotTheme->"Scientific",AxesLabel->{x,Subscript[\[Lambda], 3]},PlotLabel->"Graph of \!\(\*SubscriptBox[\(\[Lambda]\), \(3\)]\)(\!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t),0,0) versus Real \!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t)"]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
