(* ::Package:: *)

BeginPackage["SCTST2`"]
MilProb::usage="MilProb[OmegaF,xFF,forwardBarHeight,reverseBarHeight,energy] calculates the reaction probability using Miller's method (original SCTST formulation) at a paricular energy (can be above the barrier height). All inputs are in a.u."
WagProb::usage="WagProb[OmegaF,xFF,forwardBarHeight,reverseBarHeight,energy] calculates the reaction probability using Wagner's method at a paricular energy (can be above the barrier height). All inputs are in a.u."
MilBar::usage="MilBar[OmegaF,xFF,forwardBarHeight] returns the functional form of the potential barrier assumed in SCTST. All inputs are in a.u."
WagBar::usage="WagBar[OmegaF,xFF,forwardBarHeight,reverseBarHeight] returns the functional form of the potential barrier used in deep tunnelling calculations."


Begin["`Private`"]

MilProb[OmegaF_,xFF_,forBar_,revBar_,en_]:=
	Module[{tmp,deltaE},
		deltaE = forBar - en;
		tmp = OmegaF^2 + 4 xFF deltaE;
		If[And[tmp >= 0, en >= Max[forBar - revBar, 0]], (1 + Exp[2 \[Pi] (-OmegaF + Sqrt[tmp])/(2 xFF)])^-1,0]]

MilBar[OmegaF_,xFF_,forBar_]:=
	Module[{Dm,alpha},
		Dm=-OmegaF^2/4/xFF;
		alpha=Sqrt[-8 xFF];
		Function[y,forBar-Dm+4*Dm*Exp[alpha*y]/(1+Exp[alpha*y])^2]]

thetaReg[iseg_,zlo_,zhi_,eta_,rho_]:=
	Module[{c0,c1,c2,c3,c4,c5,c6,arg1l,arg2l,arg3l,thtal,arg1h,arg2h,
			arg3h,thtah},
		c0=rho^2+eta-1;
		c3=1+rho;
		c2=Sqrt[c0];
		c1=Sqrt[eta];
		If[eta>1,\[Pi](c3-c2-c1),
		c4=Sqrt[1-eta];
		c5=c3*c4;
		c6=c3-eta;
		If[iseg>=2,
			arg1l=Max[-1,Min[1,(c6*zlo-eta)/zlo/c5]];
			arg2l=Max[-1,Min[1,(c6-c0*zlo)/c5]];
			arg3l=Max[-1,Min[1,(rho*zlo-1)/(1+zlo)/c4]];
			thtal=c3*ArcSin[arg3l]+c2*ArcSin[arg2l]-c1*ArcSin[arg1l];,
			thtal=-\[Pi]/2*(c3-c2-c1);];
		If[iseg<=2,
			arg1h=Max[-1,Min[1,(c6*zhi-eta)/zhi/c5]];
			arg2h=Max[-1,Min[1,(c6-c0*zhi)/c5]];
			arg3h=Max[-1,Min[1,(rho*zhi-1)/(1+zhi)/c4]];
			thtah=c3*ArcSin[arg3h]+c2*ArcSin[arg2h]-c1*ArcSin[arg1h];,
			thtah=\[Pi]/2*(c3-c2-c1);];
		(thtah-thtal)]]

thetaIrreg[zlo_,zhi_,eta_,rho_,aa_,bb_,cc_]:=
	Module[{c1,c2,c3,arg3h,arg3l,thta,arg1,arg2,arg2h,arg2l,rhi,rlo},
		c1=Sqrt[1-eta];
		c2=1-rho^2-eta;
		c3=1+rho-eta;
		arg3h=Max[-1,Min[1,(rho*zhi-1)/(1+zhi)/c1]];
		arg3l=Max[-1,Min[1,(rho*zlo-1)/(1+zlo)/c1]];
		thta=(1+rho)*(ArcSin[arg3h]-ArcSin[arg3l]);

		rhi=Max[0,cc+bb*zhi+aa*zhi^2];
		rlo=Max[0,cc+bb*zlo+aa*zlo^2];
		arg1=(2*cc+bb*zhi+2*Sqrt[cc*rhi])/zhi;
		arg2=(2*cc+bb*zlo+2*Sqrt[cc*rlo])/zlo;
		thta-=Sqrt[-eta]*Log[arg1/arg2];

		If[-c2>=0,
			arg2h=Max[-1,Min[1,(c2*zhi+c3)/(1+rho)/c1]];
			arg2l=Max[-1,Min[1,(c2*zlo+c3)/(1+rho)/c1]];
			(thta+Sqrt[-c2]*(ArcSin[arg2h]-ArcSin[arg2l])),
			arg1=2*Sqrt[aa*rhi]+2*aa*zhi+bb;
			arg2=2*Sqrt[aa*rlo]+2*aa*zlo+bb;
			(thta+Sqrt[c2]*Log[arg1/arg2])]]

WagProb[OmegaF_,xFF_,barHeight_,revBarHeight_,energy_]:=
	Module[{rho,Db,alphab,f,ur,Rr,alphar,umr1,umr2,umr,ybr,dr,up,Rp,alphap,ump1,
			ump2,ump,ybp,dp,zbr,zrb,zbp,zpb,thtaEn,ep,eta,sqrtdelld2a,zzlo,zzhi,
			zlo,zhi,aa,bb,cc,headb,headr,headp,thta,delVf,delVr,en},
		If[barHeight>revBarHeight,
			en=energy-(barHeight-revBarHeight);
			{delVf,delVr}={revBarHeight,barHeight};,
			en=energy;
			{delVf,delVr}={barHeight,revBarHeight}];

		rho=Sqrt[delVr/delVf];
		Db=-OmegaF^2/(4*xFF)*(1-1/rho+1/rho^2);
		alphab=Abs[OmegaF]*(1+rho)/rho*Sqrt[1/2/Db];
		f=Sqrt[Db/delVf];
		ur=(Sqrt[rho^2+rho+1]-(rho-1))/Max[3*f,3];
		Rr=(210*rho^2+168*(1+f)*rho*(1-rho)*ur+140*(f-rho*(1+f)^2+f*rho^2)*ur^2-
			120*f*(1+f)*(1-rho)*ur^3+105*f^2*ur^4)/(210*rho^2+336*rho*(1-rho)*f*ur+
			140*(1-4*rho+rho^2)*f^2*ur^2-240*(1-rho)*f^3*ur^3+105*f^4*ur^4);
		alphar=alphab*f*Rr;
		umr1=(Sqrt[(1-rho)^2*(1-Rr*f)^2-4*rho*(1-Rr*f^2)*(Rr-1)]+(1-rho)*(1-Rr*f))/(2*(1-Rr*f^2));
		umr2=-(Sqrt[(1-rho)^2 (1-Rr*f)^2-4*rho*(1-Rr*f^2)*(Rr-1)]+(1-rho)*(1-Rr*f))/(2*(1-Rr*f^2));
		umr=If[Or[umr1<0,umr1>ur],umr2,umr1];
		ybr=Log[(1-umr)/(rho+umr)]/alphab;
		dr=-ybr+Log[(1-f*umr)/(rho+f*umr)]/alphar;
		up=(Sqrt[rho^2+rho+1]+rho-1)/Max[3*f,3];
		Rp=(210*rho^2-168*(1+f)*rho*(1-rho)*up+140*(f-(1+f)^2*rho+f*rho^2)up^2+120*f*(1+f)*(1-rho)*
			up^3+105*f^2*up^4)/(210*rho^2-336*rho*(1-rho)*f*up+140(1-4*rho+rho^2)f^2*up^2+240*(1-rho)*
			f^3*up^3+105*f^4*up^4);
		alphap=alphab*f*Rp;
		ump1=(Sqrt[(1-rho)^2*(1-Rp*f)^2-4*rho*(1-Rp*f^2)*(Rp-1)]-(1-rho)*(1-Rp*f))/(2*(1-Rp*f^2));
		ump2=(-Sqrt[(1-rho)^2*(1-Rp*f)^2-4*rho*(1-Rp*f^2)*(Rp-1)]-(1-rho)*(1-Rp*f))/(2*(1-Rp*f^2));
		ump=If[Or[ump1<0,ump1>up],ump2,ump1];
		ybp=Log[(-1-ump)/(-rho+ump)]/alphab;
		dp=Log[(-1-f*ump)/(-rho+f*ump)]/alphap-ybp;
		zrb=Exp[alphar*(ybr+dr)];
		zbr=Exp[alphab*ybr];
		zpb=Exp[alphap*(ybp+dp)];
		zbp=Exp[alphab*ybp];
		headb=2*Db/OmegaF*rho/(1+rho);
		headr=headb*delVf/Db/Rr;
		headp=headr*Rr/Rp;

		thtaEn=If[en<delVf,en,2*delVf-en];
		ep=thtaEn-delVf+Db;
		eta=ep/Db;
		If[thtaEn>(delVf-Db),
			sqrtdelld2a=(1+rho)*Sqrt[1-eta]/(1-rho^2-eta);
			zzlo=-(1+rho-eta)/(1-rho^2-eta)+sqrtdelld2a;
			zzhi=zzlo-2*sqrtdelld2a;
			zlo=Max[zbr,zzlo];
			zhi=Min[zbp,zzhi];
			thta=headb*thetaReg[2,zlo,zhi,eta,rho];,
			zlo=zbr;
			zhi=zbp;
			bb=2*Db*(1+rho-eta);
			aa=Db*(1-rho^2-eta);
			cc=-ep;
			thta=headb*thetaIrreg[zlo,zhi,eta,rho,aa,bb,cc];];
		If[zlo==zbr,
			eta=thtaEn/delVf;
			sqrtdelld2a=(1+rho)*Sqrt[1-eta]/(1-rho^2-eta);
			zlo=-(1+rho-eta)/(1-rho^2-eta)+sqrtdelld2a;
			zzhi=zlo-2*sqrtdelld2a;
			thta+=headr*thetaReg[1,zlo,zrb,eta,rho];
			If[zhi==zbp, thta+=headp*thetaReg[3,zpb,zzhi,eta,rho]];];
		If[en<delVf,(1+Exp[2*thta])^-1,If[(en-delVf)<delVf,1-(1+Exp[2*thta])^-1,1]]
		]

WagBar[OmegaF_,xFF_,barHeight_,revBarHeight_]:=
	Module[{rho,Db,alphab,f,ur,Rr,alphar,umr1,umr2,umr,ybr,dr,up,Rp,alphap,ump1,
			ump2,ump,ybp,dp,Vr,Vp,Vb,VWag},rho=Sqrt[revBarHeight/barHeight];
		Db=-OmegaF^2/(4*xFF)*(1-1/rho+1/rho^2);
		alphab=Abs[OmegaF]*(1+rho)/rho*Sqrt[1/2/Db];
		f=Sqrt[Db/barHeight];
		ur=(Sqrt[rho^2+rho+1]-(rho-1))/Max[3*f,3];
		Rr=(210*rho^2+168*(1+f)*rho*(1-rho)*ur+140*(f-rho*(1+f)^2+f*rho^2)*ur^2-
			120*f*(1+f)*(1-rho)*ur^3+105*f^2*ur^4)/(210*rho^2+336*rho*(1-rho)*f*ur+
			140*(1-4*rho+rho^2)*f^2*ur^2-240*(1-rho)*f^3*ur^3+105*f^4*ur^4);
		alphar=alphab*f*Rr;
		umr1=(Sqrt[(1-rho)^2*(1-Rr*f)^2-4*rho*(1-Rr*f^2)*(Rr-1)]+(1-rho)*(1-Rr*f))/(2*(1-Rr*f^2));
		umr2=-(Sqrt[(1-rho)^2 (1-Rr*f)^2-4*rho*(1-Rr*f^2)*(Rr-1)]+(1-rho)*(1-Rr*f))/(2*(1-Rr*f^2));
		umr=If[Or[umr1<0,umr1>ur],umr2,umr1];
		ybr=Log[(1-umr)/(rho+umr)]/alphab;
		dr=-ybr+Log[(1-f*umr)/(rho+f*umr)]/alphar;
		up=(Sqrt[rho^2+rho+1]+rho-1)/Max[3*f,3];
		Rp=(210*rho^2-168*(1+f)*rho*(1-rho)*up+140*(f-(1+f)^2*rho+f*rho^2)up^2+120*f*(1+f)*(1-rho)*
			up^3+105*f^2*up^4)/(210*rho^2-336*rho*(1-rho)*f*up+140(1-4*rho+rho^2)f^2*up^2+240*(1-rho)*
			f^3*up^3+105*f^4*up^4);
		alphap=alphab*f*Rp;
		ump1=(Sqrt[(1-rho)^2*(1-Rp*f)^2-4*rho*(1-Rp*f^2)*(Rp-1)]-(1-rho)*(1-Rp*f))/(2*(1-Rp*f^2));
		ump2=(-Sqrt[(1-rho)^2*(1-Rp*f)^2-4*rho*(1-Rp*f^2)*(Rp-1)]-(1-rho)*(1-Rp*f))/(2*(1-Rp*f^2));
		ump=If[Or[ump1<0,ump1>up],ump2,ump1];
		ybp=Log[(-1-ump)/(-rho+ump)]/alphab;
		dp=Log[(-1-f*ump)/(-rho+f*ump)]/alphap-ybp;
		Vr[s_]= barHeight ((1 - rho^2)*Exp[alphar*(s + dr)]/(1 + Exp[alphar*(s+dr)])+
			(1+rho)^2 *Exp[alphar*(s+dr)]/(1+Exp[alphar*(s+dr)])^2);
		Vp[s_]= barHeight *((1-rho^2)*Exp[alphap*(s+dp)]/(1+Exp[alphap*(s+dp)])+(1+rho)^2 *
			Exp[alphap*(s+dp)]/(1+Exp[alphap*(s+dp)])^2);
		Vb[s_]= barHeight-Db+Db*((1-rho^2)*Exp[alphab*s]/(1+Exp[alphab*s])+(1+rho)^2 *
			Exp[alphab*s]/(1+Exp[alphab*s])^2);
		VWag[s_]=Piecewise[{{Vr[s],s<=ybr},{Vb[s],And[s>ybr,s<ybp]},{Vp[s],s>=ybp}}];
		VWag]

End[]
EndPackage[]
