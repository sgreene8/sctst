(* ::Package:: *)

BeginPackage["JacobiPES`"]
anharm::usage"anharm[Surface,rhoSadGuess,delSadGuess,mVec,dim] calculates the frequencies and anharmonic constants (in a.u.) from an analytical surface, specified in hyperspherical coordinates. Surface is the surface function, rhoSadGuess and delSadGuess are the initial gueses for the coordinates of the saddle point on the surface. mVec is a list of the masses m1, m2, and m3, in Da. dim, an integer (1 or 2), specifies 1-D or 2-D."
surfBar::usage="surfBar[Surface,rhoMax,delProdGuess,delReacGuess,rhoSadGuess,delSadGuess,mVec] calculates the semiclassical barrier height, in a.u., for the inputted Surface. mVec, rhoSadGuess, and delSadGuess are defined as in surfAnh, rhoMax is the value of rho you want to use to find the product and reactant ZPE, and delProdGuess and delReacGuess are the initial guesses for the delta coordinates of the minima of the PES at rhoMax. In calculating the product frequency, it is assumed that the mass of the atom being abstracted is 1 Da."
MEP::usage="MEP[Surface,rhoSadGuess,delSadGuess,rhoCutoff,mVec,opt]. If opt is poten, calculates the potential (a.u.) as a function of the reaction coordinate (a.u.). Results are outputted as a 2-D table {reaction coordinate, V}. The potential is defined such that its maximum is 0. If opt is coords, outputs the hyperspherical coordinates along the MEP as {rho (a.u.), delta (radians)}, and also prints the hyperspherical coordinates of the TS. rhoCutoff is the value of rho at which the MEP should terminate. mVec is m1, m2, and m3, from the surface."


Begin["`Private`"]
Needs["PhysicalConstants`"]		

anharm[V_,rhoSadGuess_,delSadGuess_,mVec_List,dim_Integer]:=
	Module[{dVddel,dVdrho,soln,rhoSad,delSad,mu,mMat,ddr,ddR,mp,hess,rho,del,
			fMWC,tmp,modes,freqAU,ddqi,firstD,secD,thirdD,f111,f111Eval,
			f1111,f1111Eval,f2222,f1122,thirdDEval,f2222Eval,f1122Eval,anh,G0},
		dVddel[rho_,del_]=D[V[rho,del],del];
		dVdrho[rho_,del_]=D[V[rho,del],rho];
		soln=FindRoot[{dVddel[rho,del]==0,dVdrho[rho,del]==0},{rho,rhoSadGuess},{del,delSadGuess}];
		rhoSad=rho/.soln;
		delSad=del/.soln;
		mu=(mVec[[1]]*mVec[[2]]*mVec[[3]])^(1/3);
		ddr[f_]:=Sqrt[mVec[[2]]/mu]*(Sin[del]*D[f,rho]+Cos[del]/rho*D[f,del]);
		ddR[f_]:=Sqrt[mVec[[1]]/mu]*(Cos[del]*D[f,rho]-Sin[del]/rho*D[f,del]);
		mp=Convert[Dalton/ElectronMass,1];
		mMat=DiagonalMatrix[1/Sqrt[{mVec[[1]],mVec[[2]]}]];
		hess={{ddR[ddR[V[rho,del]]],ddR[ddr[V[rho,del]]]},
			{ddR[ddr[V[rho,del]]],ddr[ddr[V[rho,del]]]}};
		fMWC=mMat.hess.mMat;
		{rho,del}={rhoSad,delSad};
		{tmp,modes}=Eigensystem[fMWC];
		If[tmp[[1]]>0,
			tmp=Reverse[tmp];
			modes=Reverse[modes];];
		freqAU=Sqrt[tmp/mp];
		
		Clear[rho,del];
		ddqi[f_,i_]:={ddR[f],ddr[f]}.{modes[[i,1]]/Sqrt[mVec[[1]]],modes[[i,2]]/Sqrt[mVec[[2]]]};
		
		If[dim==1,
			f111=ddqi[ddqi[ddqi[V[rho,del],1],1],1];
			f1111=ddqi[f111,1];
			{rho,del}={rhoSad,delSad};
			f111Eval=f111/mp^(3/2);
			f1111Eval=f1111/mp^2;
			anh=1/16/freqAU[[1]]^2*(f1111Eval-5*f111Eval^2/3/freqAU[[1]]^2);
			G0=f1111Eval/64/freqAU[[1]]^2-7*f111Eval^2/576/freqAU[[1]]^4;,
			secD=modes.fMWC.Transpose[modes];
			thirdD=Table[ddqi[secD[[i,j]],k],{i,2},{j,1,i},{k,1,j}];
			f1111=ddqi[thirdD[[1,1,1]],1];
			f2222=ddqi[thirdD[[2,2,2]],2];
			f1122=ddqi[thirdD[[2,2,1]],1];
			{rho,del}={rhoSad,delSad};
			thirdDEval=thirdD/mp^(3/2);
			f1111Eval=f1111/mp^2;
			f2222Eval=f2222/mp^2;
			f1122Eval=f1122/mp^2;
			anh=Table[0,{2},{2}];
			anh[[1,1]]=1/16/freqAU[[1]]^2*(f1111Eval-5*thirdDEval[[1,1,1]]^2/3/freqAU[[1]]^2-
				thirdDEval[[2,1,1]]^2*(8*freqAU[[1]]^2-3*freqAU[[2]]^2)/freqAU[[2]]^2
				/(4*freqAU[[1]]^2-freqAU[[2]]^2));
			anh[[1,2]]=1/4/freqAU[[1]]/freqAU[[2]]*(f1122Eval-thirdDEval[[1,1,1]]*
				thirdDEval[[2,2,1]]/freqAU[[1]]^2-thirdDEval[[2,1,1]]*thirdDEval[[2,2,2]]/
				freqAU[[2]]^2+2*thirdDEval[[2,1,1]]^2*freqAU[[2]]^2/((freqAU[[1]]+
				freqAU[[2]])^2-freqAU[[1]]^2)/((freqAU[[1]]-freqAU[[2]])^2-freqAU[[1]]^2)+
				2*thirdDEval[[2,2,1]]^2*freqAU[[1]]^2/((freqAU[[1]]+freqAU[[2]])^2-
				freqAU[[2]]^2)/((freqAU[[1]]-freqAU[[2]])^2-freqAU[[2]]^2));
			anh[[2,1]]=anh[[1,2]];
			anh[[2,2]]=1/16/freqAU[[2]]^2*(f2222Eval-5*thirdDEval[[2,2,2]]^2/3/freqAU[[2]]^2-
				thirdDEval[[2,2,1]]^2*(8*freqAU[[2]]^2-3*freqAU[[1]]^2)/freqAU[[1]]^2/
				(4*freqAU[[2]]^2-freqAU[[1]]^2));
			G0=f1111Eval/64/freqAU[[1]]^2+f2222Eval/64/freqAU[[2]]^2-7*thirdDEval[[1,1,1]]^2/
				576/freqAU[[1]]^4-7/576*thirdDEval[[2,2,2]]^2/freqAU[[2]]^4+3/64*
				thirdDEval[[2,1,1]]^2/(4*freqAU[[1]]^2-freqAU[[2]]^2)/freqAU[[1]]^2+
				3/64*thirdDEval[[2,2,1]]^2/(4*freqAU[[2]]^2-freqAU[[1]]^2)/freqAU[[2]]^2;];
		{freqAU,anh,G0}]/;Length[mVec]==3&&Or[dim==1,dim==2]


surfBar[V_,rhoMax_,delProdGuess_,delReacGuess_,rhoSadGuess_,delSadGuess_,mVec_List]:=
	Module[{del,dVddel,delProd,delReac,muReac,mp,ddr,ddR,ddr1,mu,Hmass,
		rho,dVdr1dr1,dVdr2dr2,freqReac,freqProd,dVdrho,soln,rhoSad,delSad},
		dVddel[rho_,del_]=D[V[rho,del],del];
		delProd=del/.FindRoot[dVddel[rhoMax,del],{del,delProdGuess}];
		delReac=del/.FindRoot[dVddel[rhoMax,del],{del,delReacGuess}];
		mu=(mVec[[1]]*mVec[[2]]*mVec[[3]])^(1/3);
		ddr[f_]:=Sqrt[mVec[[2]]/mu]*(Sin[del]*D[f,rho]+Cos[del]/rho*D[f,del]);
		ddR[f_]:=Sqrt[mVec[[1]]/mu]*(Cos[del]*D[f,rho]-Sin[del]/rho*D[f,del]);
		mp=Convert[Dalton/ElectronMass,1];
		Hmass=ElementData["Hydrogen","AtomicWeight"];
		ddr1[f_]:=ddr[f]+mVec[[2]]/Hmass*ddR[f];
		dVdr1dr1=ddr1[ddr1[V[rho,del]]];
		dVdr2dr2=ddR[ddR[V[rho,del]]];
		
		rho=rhoMax;
		del=delProd;
		freqProd=Sqrt[dVdr1dr1/mVec[[2]]/mp];
		
		del=delReac;
		muReac=mVec[[1]]*Hmass^2/(mVec[[1]]*mVec[[2]]+Hmass^2)*mp;
		freqReac=Sqrt[dVdr2dr2/muReac];

		dVdrho[rho_,del_]=D[V[rho,del],rho];
		soln=FindRoot[{dVddel[rho,del]==0,dVdrho[rho,del]==0},{rho,rhoSadGuess},{del,delSadGuess}];
		rhoSad=rho/.soln;
		delSad=del/.soln;
		
		{V[rhoSad,delSad]-V[rhoMax,delReac]-freqReac/2,V[rhoSad,delSad]-V[rhoMax,delProd]-freqProd/2}]/;Length[mVec]==3

MEP[V_,rhoSadGuess_,delSadGuess_,rhoCutoff_,mVec_List,opt_String]:=
	Module[{dVddel,dVdrho,soln,rhoSad,delSad,mu,mMat,ddr,ddR,hess,fMWC,direction,dVdr,dVdR,i,
			rho,del,tmp,modes,rxnMode,delq=0.001,VSDProd,VSDReac,Ri,ri,coordsProd,coordsReac,
			grad,maxPos,mp},
		dVddel[rho_,del_]=D[V[rho,del],del];
		dVdrho[rho_,del_]=D[V[rho,del],rho];
		soln=FindRoot[{dVddel[rho,del]==0,dVdrho[rho,del]==0},{rho,rhoSadGuess},{del,delSadGuess}];
		rhoSad=rho/.soln;
		delSad=del/.soln;
		mu=(mVec[[1]]*mVec[[2]]*mVec[[3]])^(1/3);
		ddr[f_]:=Sqrt[mVec[[2]]/mu]*(Sin[del]*D[f,rho]+Cos[del]/rho*D[f,del]);
		ddR[f_]:=Sqrt[mVec[[1]]/mu]*(Cos[del]*D[f,rho]-Sin[del]/rho*D[f,del]);
		dVdr=ddr[V[rho,del]];
		dVdR=ddR[V[rho,del]];
		mMat=DiagonalMatrix[1/Sqrt[{mVec[[1]],mVec[[2]]}]];
		hess={{ddR[dVdR],ddR[dVdr]},
			{ddR[dVdr],ddr[dVdr]}};
		fMWC=mMat.hess.mMat;
		{rho,del}={rhoSad,delSad};
		{tmp,modes}=Eigensystem[fMWC];
		If[tmp[[1]]>0,modes=Reverse[modes];];
		
		VSDProd=Table[0,{5000}];
		coordsProd=Table[{0,0},{5000}];
		coordsProd[[1]]={rho,del};
		VSDProd[[1]]=V[rhoSad,delSad];
		{Ri,ri}=rho*Sqrt[mu]*{Cos[del]/Sqrt[mVec[[1]]],Sin[del]/Sqrt[mVec[[2]]]};
		direction=If[modes[[1,2]]<0,1,-1];
		{Ri,ri}+=delq*direction*modes[[1]]*mVec[[{1,2}]]^-.5;
		For[i=2,And[i<=5000,rho<=rhoCutoff],i++,
			rho=Sqrt[mVec[[1]]/mu Ri^2+mVec[[2]]/mu ri^2];
			del=ArcTan[Sqrt[mVec[[2]]/mVec[[1]]] ri/Ri];
			coordsProd[[i]]={rho,del};
			grad=mVec[[{1,2}]]^-.5*{dVdR,dVdr};
			VSDProd[[i]]=V[rho,del];
			{Ri,ri}-=delq*mVec[[{1,2}]]^-.5*grad/Norm[grad];];
		(*Get rid of extra zeros*)
		VSDProd=Delete[VSDProd,Position[VSDProd,0]];
		coordsProd=Delete[coordsProd,Position[coordsProd,{0,0}]];

		VSDReac=Table[0,{5000}];
		coordsReac=Table[{0,0},{5000}];
		{rho,del}={rhoSad,delSad};
		coordsReac[[1]]={rho,del};
		VSDReac[[1]]=V[rhoSad,delSad];
		{Ri,ri}=rho*Sqrt[mu]*{Cos[del]/Sqrt[mVec[[1]]],Sin[del]/Sqrt[mVec[[2]]]};
		{Ri,ri}-=delq*direction*modes[[1]]*mVec[[{1,2}]]^-.5;
		For[i=2,And[i<=5000,rho<=rhoCutoff],i++,
			rho=Sqrt[mVec[[1]]/mu Ri^2+mVec[[2]]/mu ri^2];
			del=ArcTan[Sqrt[mVec[[2]]/mVec[[1]]] ri/Ri];
			coordsReac[[i]]={rho,del};
			grad=mVec[[{1,2}]]^-.5*{dVdR,dVdr};
			VSDReac[[i]]=V[rho,del];
			{Ri,ri}-=delq*mVec[[{1,2}]]^-.5*grad/Norm[grad];];
		VSDReac=Delete[VSDReac,Position[VSDReac,0]];
		coordsReac=Delete[coordsReac,Position[coordsReac,{0,0}]];

		If[opt=="poten",
			tmp=Join[Reverse[VSDReac],Delete[VSDProd,1]];
			maxPos=Position[tmp,Max[tmp]][[1,1]];
			mp=Convert[Dalton/ElectronMass,1];
			Table[{delq*Sqrt[mp]*(i-maxPos),tmp[[i]]-V[rhoSad,delSad]},{i,Length[tmp]}],
			
			If[opt=="coords",Join[Reverse[coordsReac],Delete[coordsProd,1]],0]]
		]/;Length[mVec]==3

End[]
EndPackage[]
