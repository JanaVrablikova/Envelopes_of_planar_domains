(* ::Package:: *)

BeginPackage["Functions`"];

MInnProd::usage="MInnProd[v1_,v2] returns the Minkowski inner product of two vectors in R^(2,1)."
MSignatureMat::usage="returns the signature matrix of R^(2,1)."

BCurve::usage="BCurve[pts_,var_] returns a Bezier curve with control points -pts- in -var-."
BSurf::usage="BCurve[pts_,var1_,var2_] returns a Bezier surface with control points -pts- in -var1- and -var2-."

Env::usage="Env[c_,var_,tol_] returns the envelope of the Minkowski curve -c- in -var-. The curve is expected to be parametrized over [0,1]. The envelope is parametrized in -var- over [0,1]. If length(-c-) < -tol-, it returns the cyclographic image of -c(0)-)." 

Newton::usage="Newton[f_,var_,var0_,tol_,maxit_] minimizes -f-(-var-) starting from -var0- with tolerance -tol- and at most -maxit- iterations."

MinkArcParam::usage="ArcParams[pts_,var_,tol_] returns a Minkowski arc parametrized in -var- over [0,1]. The arc interpolates the Minkowski points -pts-."

BezierArc::usage="BezierArc[pts_,sign_,var_] returns circular arc with three control points -pts- parametrized in -var- over [0,1]. If -sign-=1, the arc lies in the convex hull ot the points, if -sign-=-1, it returns its complement."

CPtsFromCentre::usage="[p1_,p3_,s_,tg1_] returns the middle control point and a -sign- of a Bezier arc interpolating points -p1- and -p3- with -s- being the centre of the arc and -tg1- being the tangent vector of the arc at -p1-. sign in 1 if the arc lies in the convex hull of the control points, -1 otherwise."

ArcFrom3Pts::usage="[p1_,p2_,p3_,var_] returns the circular arc interpolating the three points -p1-,-p2-,p3- parametrized in -var- over [0,1]."
 
MinkBiarc::usage="MinkBiarc[p1_,tg1_,p2_,tg2_,var_,tol_] returns a Minkowski biarc interpolating the data (-p1-,-tg1-), (-p2-,-tg2-), parametrized over [0,1] in -var-."

BiArc::usage="BiArc[p1_,tg1_,p2_,tg2_,var_,tol_] returns a biarc interpolating the data (-p1-,-tg1-), (-p2-,-tg2-), parametrized over [0,1] in -var-."

EndCap::usage="EndCap[c_,var_,tol_] returns arcs of the cyclographic images of -c-(0) and -c-(1), such that the envelope of -c- together with the arcs form the boundary of the cyclographic image of the curve -c-. -var- is the parameter of the curve."


Begin["`Private`"];

MSignatureMat={{1,0,0},{0,1,0},{0,0,-1}};

MInnProd[pt1_,pt2_]:=(pt1 . MSignatureMat) . pt2

mDist2[pt1_,pt2_]:=MInnProd[pt1-pt2,pt1-pt2]


bBasis[var_,n_]:=Table[Binomial[n,i]*var^i*(1-var)^(n-i),{i,0,n}]
BCurve[pts_,t_]:=Sum[pts[[i]]*bBasis[t,Length@pts-1][[i]],{i,1,Length@pts}]
BSurf[pts_,u_,t_]:=Sum[Sum[pts[[i,j]]*bBasis[u,Length@pts-1][[i]]*bBasis[t,Length@pts[[i]]-1][[j]],{j,1,Length@pts[[i]]}],{i,1,Length@pts}]


Env[c_,var_,tol_]:=If[NIntegrate[Total[D[c,var]^2],{var,0,1},AccuracyGoal->6]<tol,
	Return[{(c[[1;;2]]/.{var->0})+(c[[3]]/.{var->0})*{Cos[Pi*var],Sin[Pi*var]},(c[[1;;2]]/.{var->0})-(c[[3]]/.{var->0})*{Cos[Pi*var],Sin[Pi*var]}}]
,
	Return[Table[c[[1;;2]]-c[[3]]*(D[c[[3]],var]*D[c[[1;;2]],var]-(-1)^i*Sqrt[MInnProd[D[c,var],D[c,var]]]*{D[c[[2]],var],-D[c[[1]],var]} )/({D[c[[1]],var],D[c[[2]],var]}.{D[c[[1]],var],D[c[[2]],var]}),{i,1,2}]]
];


Newton[f_,var_,var0_,tol_,maxit_]:=
Module[{x0,fx0,dfx0,x1,i,xerr},(
	x0=var0;
	For[i=1,i<=maxit,i++,(
		fx0=f/.{var->x0};
		dfx0=((f/.{var->x0+tol/2})-(f/.{var->x0-tol/2}))/(tol);
		If[Or[Im@dfx0>0,Im@fx0>0],
			Return[{0,0,0},Module]];
		If[Abs@dfx0<tol,
			Return[{0,0,0},Module]];
		x1=x0-fx0/dfx0;
		If[Abs[x1-x0]<tol,
			Return[{1,x1,f/.{var->x1}},Module]];
		xerr=x1-x0;
		x0=x1;
		(*Print[x0];*)
	)];
	(*Print["did not converge in "<>ToString[maxit]<>" itertions"<>", error = "<>ToString[Norm[xerr]]];*)
Return[{0,0,0},Module];)]


(*for three points -pts- such that |p2-p1| = |p2-p3| it reruns three wegihts (1,sign*w,1) for a rational bezier arc with these control points*)
weights[pts_,sign_]:={1,sign*Norm[(pts[[1]]-pts[[3]])/2]/Norm[pts[[1]]-pts[[2]]],1};

BezierArc[pts_,sign_,var_]:=Module[{a},
(
    If[Length@pts == 2,
    (*line segment*)
		a=pts[[1]]+var*(pts[[2]]-pts[[1]]);,
    (*else, arc*)
		a=Sum[pts[[i]]*bBasis[var,2][[i]]*weights[pts,sign][[i]],{i,1,3}]/Sum[bBasis[var,2][[i]]*weights[pts,sign][[i]],{i,1,3}]
    ];
    Return[a];
)]

CPtsFromCentre[p1_,p3_,s_,tg1_]:=Module[{r,M,sinA,p2,sign,ll},
(
	r=Norm[s-p1];
	M=(p1+p3)/2;
	sign=Sign@Dot[tg1,p3-p1];
	sinA=Norm[s-M]/r;
	p2=M+(r*(1-sinA^2)/sinA)*Normalize[M-s];
	Return[{p2,sign},Module]
)]

ArcFrom3Pts[pt1_,pt2_,pt3_,var_]:=Module[{centre,l1,s1,l2,s2,sol,cp2,sign},
(
	l1=(pt1+pt2)/2+s1*{{0,-1},{1,0}} . (pt2-pt1);
	l2=(pt2+pt3)/2+s2*{{0,-1},{1,0}} . (pt3-pt2);
	sol=NSolve[l1==l2,{s1,s2}];
	If[Length@sol==0,
		Return[(1-var)*pt1+var*pt3,Module]
		,
		(
		centre=l1/.sol[[1]];
		{cp2,sign}=CPtsFromCentre[pt1,pt3,centre,-Sign[Det[{pt2-pt1,centre-pt1}]]*{{0,-1},{1,0}} . (centre-pt1)];
		Return[BezierArc[{pt1,cp2,pt3},sign,var],Module]
		)]
)]



bisector[p1_,p2_,var_]:=Return[(p1+p2)/2+var*{{0,-1},{1,0}} . (p2-p1)]

BiArc[p1_,tg1_,p2_,tg2_,var_,tol_]:=Module[
{biarcs,cPts,s1,s2,l1,l2,Cj,J,rj,centre,radii,a1,a2,tg1n,tg2n,biarc},
(
biarcs={0,0};
cPts={0,0};
centre={0,0};
radii={0,0};
tg1n=Normalize[tg1];
tg2n=Normalize[tg2];
If[Abs@Det[{(p2-p1),(p2+tg2n-p1-tg1n)}]<tol,(
		J=(p1+p2)/2;
		cPts={{p1,J},{J,p2}};
		biarc={BezierArc[cPts[[1]],1,var],BezierArc[cPts[[2]],1,var]};
		centre={0,0};
		radii={0,0};
		Return[{cPts,biarc,centre,radii},Module]
)];
l1=bisector[p1,p2,s1];
l2=bisector[p1+tg1n,p2+tg2n,s2];
Cj=l1/.NSolve[l1==l2,{s1,s2}][[1]];
rj=Norm[p1-Cj];
a1=BezierArc[{p1,CPtsFromCentre[p1,p2,Cj,tg1n][[1]],p2},1,var];
a2=BezierArc[{p1,CPtsFromCentre[p1,p2,Cj,tg1n][[1]],p2},-1,var];
If[NIntegrate[Total[D[a1,var]^2],{var,0,1}]<NIntegrate[Total[D[a2,var]^2],{var,0,1}],
	J=a1/.{var->1/2},
	J=a2/.{var->1/2}
];
J=a1/.{var->1/2};
l1=bisector[p1,J,s1];
centre[[1]]=l1/.Solve[l1==p1+s2*({{0,-1},{1,0}} . tg1n),{s1,s2}][[1]];
radii[[1]]=Norm[centre[[1]]-p1];
cPts[[1]]={p1,CPtsFromCentre[p1,J,centre[[1]],tg1n][[1]],J};
biarcs[[1]]=BezierArc[cPts[[1]],1,var];
l2=bisector[J,p2,s1];
centre[[2]]=l2/.Solve[l2==p2+s2*({{0,-1},{1,0}} . tg2n),{s1,s2}][[1]];
radii[[2]]=Norm[centre[[2]]-p2];
cPts[[2]]={J,CPtsFromCentre[J,p2,centre[[2]],{{0,-1},{1,0}} . (Cj-centre[[2]])][[1]],p2};
biarcs[[2]]=BezierArc[cPts[[2]],1,var];
Return[{cPts,biarcs,centre,radii,Cj,rj,a1}])]


MinkArcParam[pts_,u_,tol_]:=Module[{pts1,bases,coefs,alphaN,alphaD,betaN,betaD,arc,x1,x2,y1,y2,r1,r2},
(
pts1=Table[Join[{1},pts[[i]]-pts[[1]]],{i,1,Length@pts}];
bases={(u-1)*(u-1/2),u*(u-1),u*(u-1/2)};
betaN=(r1^2-x1^2-y1^2)/.{x1->pts1[[2,2]],y1->pts1[[2,3]],r1->pts1[[2,4]],x2->pts1[[3,2]],y2->pts1[[3,3]],r2->pts1[[3,4]]};
betaD=(r1^2-2 r1 r2+r2^2-x1^2+2 x1 x2-x2^2-y1^2+2 y1 y2-y2^2)/.{x1->pts1[[2,2]],y1->pts1[[2,3]],r1->pts1[[2,4]],x2->pts1[[3,2]],y2->pts1[[3,3]],r2->pts1[[3,4]]};
alphaN=(-r2^2+x2^2+y2^2)/.{x1->pts1[[2,2]],y1->pts1[[2,3]],r1->pts1[[2,4]],x2->pts1[[3,2]],y2->pts1[[3,3]],r2->pts1[[3,4]]};
alphaD=(-2 (-r1^2+2 r1 r2-r2^2+x1^2-2 x1 x2+x2^2+y1^2-2 y1 y2+y2^2))/.{x1->pts1[[2,2]],y1->pts1[[2,3]],r1->pts1[[2,4]],x2->pts1[[3,2]],y2->pts1[[3,3]],r2->pts1[[3,4]]};
	If[Abs[betaD]<=tol||Abs[alphaD]<=tol,Return[pts[[1]],Module]];
coefs={1,alphaN/alphaD,betaN/betaD};
arc=Sum[pts1[[i]]*bases[[i]]*coefs[[i]],{i,1,3}];
arc=arc[[2;;4]]/arc[[1]];
Return[Simplify[arc+pts[[1]]]]
)]


minkWeights[pts_]:={1,Sqrt[MInnProd[(pts[[1]]-pts[[3]])/2,(pts[[1]]-pts[[3]])/2]]/Sqrt[MInnProd[pts[[1]]-pts[[2]],pts[[1]]-pts[[2]]]],1};
MbezierArc[pts_,var_]:=Module[{a},
(
    If[Length@pts == 2,
    (*line segment*)
		a=pts[[1]]+var*(pts[[2]]-pts[[1]]);,
    (*else, arc*)
		a=Sum[pts[[i]]*bBasis[var,2][[i]]*minkWeights[pts][[i]],{i,1,3}]/Sum[bBasis[var,2][[i]]*minkWeights[pts][[i]],{i,1,3}]
    ];
    Return[a];
)]

MinkBiarc[p1_,tg1_,p2_,tg2_,var_,tol_]:=Module[{b1,b2,l1,l2,eq,\[Lambda]1,\[Lambda]2,J,sol,ll1,ba1,ba2},(
b1=p1+tg1*l1;
b2=p2-tg2*l2;
eq=Simplify[mDist2[p1,p2]-2l1*MInnProd[p2-p1,tg1]-2l2*MInnProd[p2-p1,tg2]+2*l1*l2*(MInnProd[tg1,tg2]-1)];
\[Lambda]2=l2/.Solve[eq==0,l2][[1]]/.{l1->\[Lambda]1};
b1=b1/.{l1->\[Lambda]1};
b2=b2/.{l2->\[Lambda]2};
J=Simplify[(\[Lambda]2*b1+\[Lambda]1*b2)/(\[Lambda]1+\[Lambda]2)];
sol=NSolve[\[Lambda]2==\[Lambda]1,\[Lambda]1];
If[Total[Abs[Im[\[Lambda]1/.sol]]]>0,Return[{0,0},Module]];
\[Lambda]1=Sort[\[Lambda]1/.sol,Abs[#1]<Abs[#2]&][[1]];
ba1=MbezierArc[{p1,b1,J},var];
ba2=MbezierArc[{J,b2,p2},var];
Return[{ba1,ba2},Module];
)]


envPts[c_,var_,pars_,tol_]:=Table[Env[c,var,tol]/.{var->pars[[i]]},{i,1,Length@pars}]

EndCap[c_,var_,tol_]:=Module[{epts,cps,ecaps,endcaps,p1,p2,tg,tgl,sol,s,z,r,arc},(
	epts=envPts[c,var,{0,1},tol];
	endcaps={0,0};
	Table[
		If[Norm[epts[[i,1]]-epts[[i,2]]]<tol,
			endcaps[[i]]={1,1}
			,
			If[Abs@Det[Join[{epts[[i,1]],(c/.{var->(i-1)})[[1;;2]],epts[[i,2]]},{{1},{1},{1}},2]]<tol, (*half-circle*)
			(
				p1=epts[[i,1]];
				z=(c/.{var->(i-1)})[[1;;2]];
				p2=epts[[i,2]];
				r=(c/.{var->(i-1)})[[3]];
				arc=z+r*{Cos[var],Sin[var]}/.{var->var+ArcCos[Norm[p1[[1]]-z[[1]]]/Norm[p1-z]]};
				ecaps={arc/.{var->var*Pi}, arc/.{var->-var*Pi}};
				tg=Normalize[D[c[[1;;2]],var]/.{var->i-1}];
				If[Norm[(z+tg)-(ecaps[[1]]/.{var->1/2})]>Norm[(z+tg)-(ecaps[[2]]/.{var->1/2})],
					(
					If[i==1,
						endcaps[[i]]=ecaps[[1]],
						endcaps[[i]]=ecaps[[2]]
					]),
					(
					If[i==1,
						endcaps[[i]]=ecaps[[2]],
						endcaps[[i]]=ecaps[[1]]
					])
				];
			),
			(
				cps=CPtsFromCentre[epts[[i,1]],epts[[i,2]],(c/.{var->(i-1)})[[1;;2]],{{0,-1},{1,0}}.Normalize[epts[[i,1]]-(c/.{var->(i-1)})[[1;;2]]]];
				ecaps={BezierArc[{epts[[i,1]],cps[[1]],epts[[i,2]]},1,var],BezierArc[{epts[[i,1]],cps[[1]],epts[[i,2]]},-1,var]};
				tg=D[c[[1;;2]],var]/.{var->i-1};
				tgl=(c[[1;;2]]/.{var->i-1})+s*tg;
				sol={NSolve[tgl==ecaps[[1]]&&0<=var<=1&&0<=s,{var,s},Reals],NSolve[tgl==ecaps[[2]]&&0<=var<=1&&0<=s,{var,s},Reals]};
				If[Length[sol[[1]]]==0,
					If[i==1,
						endcaps[[i]]=ecaps[[1]],
						endcaps[[i]]=ecaps[[2]]
					],
					If[i==1,
						endcaps[[i]]=ecaps[[2]],
						endcaps[[i]]=ecaps[[1]]
					]
				]
			)];
		]
	,{i,1,2}];
	Return[endcaps,Module];
)]


End[];(*`Private`*)
EndPackage[];

