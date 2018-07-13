(* ::Package:: *)

(* ::Input:: *)
(*trueVanilla = {r->0.2,k->10,crit->0.1,ns->100,no->100,ps->5000,po->5000,qs->0.007,qo->0.007,ws->1,wo->1,ms->1,mo->1,ehats->0.6,ehato->0.6,gammas->0.1,gammao->0.1,betas->0.2,betao->0.2,a->-0.3,b->0.1102};*)
(**)
(*(*Graph under 0 < R_S < ehats constraint.*) *)
(*(*Pre-analysis*)*)
(*(*Derivative of profit of defectors will always be 21*)*)
(*Dpids[x_]=ps*qs*ehats/.trueVanilla;*)
(*(*Derivative of profit of cooperators will always be 1.7938000000000003\[VeryThinSpace]+0.1102 x*)*)
(*Dpics[x_]=(ps*qs*ehats-ws-ms)b+(a+x*b)/.trueVanilla;*)
(*(*ecs[x] will always be R_S*)*)
(*ecs[x_]=a+x*b/.trueVanilla;*)
(*(*Effort of defector SON*)*)
(*eds[x_]=ehats;*)
(*(*Effort of defector OMNRF*)*)
(*edo[x_]=ehato*)
(*(*Profit of Defectors SON is -0.6+21.x*)*)
(*pids[x_]=ps*qs*ehats*x-ws*ehats/.trueVanilla;*)
(*(*Profit of Coop SON is 1.2-0.2204x+35.(-0.3+0.1102 x)x*)*)
(*pics[x_]=ps*qs*(a+b*x)*x-ws*(a+x*b)+ms*(ehats-(a+x*b))/.trueVanilla;*)
(*(*CS2*)*)
(*CS2=ns*((-betas/gammas)*(1-(pics[x]/pids[x]))+1)/.trueVanilla;*)
(*Plot[CS2,{x,(-a/b),((ehats-a)/b)}/.trueVanilla]*)
(**)
(*(*EV2*)*)
(*-1/(2 k ns pids[x]^2) (-2 betas CS2 gammas k pics[x] pids[x]+(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x]^2+\[Sqrt](pids[x]^2 ((-2 betas CS2 gammas k pics[x]+(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x])^2+4 k ns (betas CS2 k ns qs x Dpids[x] (ecs[x]-eds[x]) pics[x]+pids[x] (betas CS2 k ns qs x Dpics[x] (-ecs[x]+eds[x])+gammas (k (CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) r x+3 r x^2+k no qo edo[x]+k qs (CS2 qs ecs[x]+(-CS2+ns) eds[x])) (2 betas CS2 pics[x]+(-2 betas CS2+ns) pids[x]))))))/.trueVanilla;*)
(**)
(*Plot3D[-(1/(2 k ns pids[x]^2))(-2 betas CS2 gammas k pics[x] pids[x]+(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x]^2+\[Sqrt](pids[x]^2 ((-2 betas CS2 gammas k pics[x]+(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x])^2+4 k ns (betas CS2 k ns qs x Dpids[x] (ecs[x]-eds[x]) pics[x]+pids[x] (betas CS2 k ns qs x Dpics[x] (-ecs[x]+eds[x])+gammas (k (CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) r x+3 r x^2+k no qo edo[x]+k qs (CS2 qs ecs[x]+(-CS2+ns) eds[x])) (2 betas CS2 pics[x]+(-2 betas CS2+ns) pids[x]))))))/.trueVanilla,{x,(-a/b),((ehats-a)/b)}/.trueVanilla,{eco,0,0.6},AxesLabel->{Subscript[SuperStar[x], 2][t],Subscript[(e^C), o][t],Subscript[\[Lambda], 1]},PlotLabel->"Graph of \!\(\*SubscriptBox[\(\[Lambda]\), \(1\)]\)(\!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t),\!\(\*SubscriptBox[SuperscriptBox[\(C\), \(*\)], SubscriptBox[\(s\), \(2\)]]\)(t),0) versus Real \!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t) and OMNRF Cooperator Effort from 0 to 0.6"]*)
(**)
(*(*EV3*)*)
(*1/(2 k ns pids[x]^2) (2 betas CS2 gammas k pics[x] pids[x]-(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x]^2+\[Sqrt](pids[x]^2 ((-2 betas CS2 gammas k pics[x]+(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x])^2+4 k ns (betas CS2 k ns qs x Dpids[x] (ecs[x]-eds[x]) pics[x]+pids[x] (betas CS2 k ns qs x Dpics[x] (-ecs[x]+eds[x])+gammas (k (CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) r x+3 r x^2+k no qo edo[x]+k qs (CS2 qs ecs[x]+(-CS2+ns) eds[x])) (2 betas CS2 pics[x]+(-2 betas CS2+ns) pids[x]))))))/.trueVanilla;*)
(*Plot3D[1/(2 k ns pids[x]^2) (2 betas CS2 gammas k pics[x] pids[x]-(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x]^2+\[Sqrt](pids[x]^2 ((-2 betas CS2 gammas k pics[x]+(2 betas CS2 gammas k+k ns (-gammas+CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) ns r x+3 ns r x^2+k ns (CS2 qs^2 ecs[x]+no qo edo[x]+(-CS2+ns) qs eds[x])) pids[x])^2+4 k ns (betas CS2 k ns qs x Dpids[x] (ecs[x]-eds[x]) pics[x]+pids[x] (betas CS2 k ns qs x Dpics[x] (-ecs[x]+eds[x])+gammas (k (CS2 eco qo^2+b CS2 qs+crit r)-2 (crit+k) r x+3 r x^2+k no qo edo[x]+k qs (CS2 qs ecs[x]+(-CS2+ns) eds[x])) (2 betas CS2 pics[x]+(-2 betas CS2+ns) pids[x]))))))/.trueVanilla,{x,(-a/b),((ehats-a)/b)}/.trueVanilla,{eco,0,0.6},AxesLabel->{Subscript[SuperStar[x], 2][t],Subscript[(e^C), o][t],Subscript[\[Lambda], 2]},PlotLabel->"Graph of \!\(\*SubscriptBox[\(\[Lambda]\), \(2\)]\)(\!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t),\!\(\*SubscriptBox[SuperscriptBox[\(C\), \(*\)], SubscriptBox[\(s\), \(2\)]]\)(t),0) versus Real \!\(\*SubscriptBox[SuperscriptBox[\(x\), \(*\)], \(2\)]\)(t) and OMNRF Cooperator Effort from 0 to 0.6"]*)


(* ::Input:: *)
(*Show[%86,Boxed->False]*)


(* ::Input:: *)
(*Show[%68,Boxed->False]*)


(* ::Input:: *)
(*Show[%71]*)
