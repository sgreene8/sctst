(* ::Package:: *)

BeginPackage["AbIniAnh`"]
adDiff::usage="adDiff[numDerivs,stepFactor] reads in a list of numerical derivatives. Each derivative in the list is assumed to be calculated with a step size greater than the previous value by a factor stepFactor." 
vibAnalysis::usage="vibAnalysis[paths,masses,opt,atomNames,nHess,nQ,dQ,outDir,memory] reads in a list of paths to hessian output files calculated with increasing step sizes (by a factor of 2). If opt is 'AD results,' performs AD 
on all hessian elements and outputs the resulting tensor. If opt is 'freqs,' outputs the frequencies from the most accurate converged hessian, using the list masses of atomic masses. If opt is 'Input files,' calculates the displaced geometries from 
the converged Hessian and creates input files for MOLPRO in the directory outDir. atomNames is a list of the names of all of the atoms (C, H, etc.), nHess is the number of hessian calculations to perform at each displaced geom (smallest step size is 0.01), 
nQ is the number of deltaQ steps, dQ sets the smallest deltaQ step size in reduced normal coordinates, memory specifies the amount of memory to use for each hessian calculation. If opt is '1-D en,' calculates the displaced geometries along the reaction
mode and creates input files for a number of deltaQ steps specified by nQ with a step size of dQ in reduced normal coordinates."
calcAnh::usage="calcAnh[hessPaths,masses,anhDir,nHess,nQ,dQ,nModes] calculates the anharmonic constants from ab initio hessian calculations. hessPaths is a list of hessian output files for the TS (calculated with step sizes increasing by a factor of 2, 
masses is the list of atomic masses, anhDir is the path of the directory where the displaced geometry hessian output files are located, nHess is the number of hessian matrices calculated at each displaced geometry, nQ is the number of deltaQ step sizes
used, dQ is the smallest deltaQ step size in reduced normal coordinates, and nModes specifies the number of modes to be treated."
anh1D::usage="anh1D[energies,dQ,freq] calculates the 1-D anharmonic constant and G0 from single point energies along the reaction mode. The List energies should have an odd number of elements, with the TS energy as the middle element, and the displacements 
from the TS increasing by factors of 2."


Begin["`Private`"]
Needs["readInTxt`"]
Needs["PhysicalConstants`"]

adDiff[numDerivs_List,stepFactor_]:=
	Module[{outMat,fac,n=Length[numDerivs],k,m},
		outMat=Table[0,{n},{n}];
		For[k=1,k<=n,k++,
			outMat[[1,k]]=Reverse[numDerivs][[k]];
			fac=1;
			For[m=2,m<=k,m++,
				fac*=stepFactor^2;
				outMat[[m,k]]=(outMat[[m-1,k]]*fac-outMat[[m-1,k-1]])/(fac-1);]];
		outMat]

printCoordsHess[coords_List,names_List,fName_String,stepSize_,memory_Integer]:=
	Module[{file,i,numbers},
		file=OpenWrite[fName,FormatType->StandardForm];
		WriteString[file,"***,",fName,"\n\nmemory,"<>IntegerString[memory]<>",m\n\nnosym\nangstrom\ngeometry={\n"];
		For[i=1,i<=Length[coords],i++,
			numbers=Table[NumberForm[coords[[i,k]]*Convert[BohrRadius/Angstrom,1],{11,10},
				ExponentFunction->(Null &)],{k,3}];
			Write[file," ",names[[i]],"   ",numbers[[1]],"   ",numbers[[2]],"   ",numbers[[3]]]];
		WriteString[file,"}\n\nbasis=vtz\nuhf,maxit=400\nump2,maxit=400\n"];
		WriteString[file,"{frequencies,step=0."<>IntegerString[IntegerPart[100*stepSize],10,2]
			<>",PRINT=1}\n"];
		WriteString[file,"\n---\n"];
		Close[file];]/;Length[coords]==Length[names]

printCoordsEn[coords_List,names_List,fName_String,memory_Integer]:=
	Module[{file,i,numbers},
		file=OpenWrite[fName,FormatType->StandardForm];
		WriteString[file,"***,",fName,"\n\nmemory,"<>IntegerString[memory]<>",m\n\nnosym\nangstrom\ngeometry={\n"];
		For[i=1,i<=Length[coords],i++,
			numbers=Table[NumberForm[coords[[i,k]]*Convert[BohrRadius/Angstrom,1],{11,10},
				ExponentFunction->(Null &)],{k,3}];
			Write[file," ",names[[i]],"   ",numbers[[1]],"   ",numbers[[2]],"   ",numbers[[3]]]];
		WriteString[file,"}\n\nbasis=vtz\nuhf,maxit=400\nump2,maxit=400\n\n---\n"];
		Close[file];]/;Length[coords]==Length[names]

vibAnalysis[fNames_List,masses_List,opt_String,atoms_List,nHess_Integer,nQ_Integer,dQ_,outDir_String,memory_Integer]:=
	Module[{mMat,mVec,hessData,hessCalc,n=Length[masses],s=Length[fNames],hessConv,
			mwh,vals,modes,freqs,mp,freqConv,tmp,tsGeom,Qsteps,a,b,c,cartVec,pgeom,mgeom,dx},
		hessData=Table[readInTxt[fNames[[i]],n,0][[2]],{i,s}];
		hessCalc=Table[If[i>=j,adDiff[hessData[[All,i,j]],2],
			Table[0,{s},{s}]],{i,3*n},{j,3*n}];
		If[opt=="AD results",hessCalc,
		hessConv=hessCalc[[All,All,s,s]];
		hessConv=hessConv+Transpose[hessConv]-DiagonalMatrix[Diagonal[hessConv]];
		mVec=Table[masses[[Floor[(i-1)/3]+1]],{i,3*n}];
		mMat=DiagonalMatrix[1/Sqrt[mVec]];
		mwh=mMat.hessConv.mMat;
		{vals,modes}=Eigensystem[mwh,3*n-6];
		mp=Convert[Dalton/ElectronMass,1];
		freqConv=Convert[Centi*Meter/(2*\[Pi]*SpeedOfLight)*Sqrt[2*Rydberg/BohrRadius^2/Dalton],1];
		freqs=Sqrt[vals/mp];
		tmp=Position[Re[freqs],0.][[1]];
		freqs=Join[freqs[[tmp]],Delete[freqs,tmp]];
		modes=Join[modes[[tmp]],Delete[modes,tmp]];
		If[opt=="freqs",freqs*freqConv*Sqrt[mp],
		tsGeom=readInTxt[fNames[[1]],n,0][[1]];
		Qsteps=dQ/Sqrt[Abs[freqs]*mp];
		If[opt=="1-D en",
		cartVec=Table[(modes[[1]]/Sqrt[mVec])[[3(i-1)+j]],{i,n},{j,3}];
		For[b=1,b<=nQ,b++,
			tmp=Qsteps[[1]]*2^(b-1);
			pgeom=tsGeom+tmp*cartVec;
			mgeom=tsGeom-tmp*cartVec;
			printCoordsEn[pgeom,atoms,outDir<>"/p"<>IntegerString[b]<>".txt",memory];
			printCoordsEn[mgeom,atoms,outDir<>"/m"<>IntegerString[b]<>".txt",memory];];,
		If[opt=="Input files",
		For[a=1,a<=(3*n-6),a++,
			cartVec=Table[(modes[[a]]/Sqrt[mVec])[[3(i-1)+j]],{i,n},{j,3}];
			For[b=1,b<=nQ,b++,
				tmp=Qsteps[[a]]*2^(b-1);
				pgeom=tsGeom+tmp*cartVec;
				mgeom=tsGeom-tmp*cartVec;
				For[c=1,c<=nHess,c++,
					dx=0.01*2^(c-1);
					printCoordsHess[pgeom,atoms,outDir<>"/"<>IntegerString[b]<>"p"<>IntegerString[a]<>
						IntegerString[IntegerPart[100*dx],10,2]<>".txt",dx,memory];
					printCoordsHess[mgeom,atoms,outDir<>"/"<>IntegerString[b]<>"m"<>IntegerString[a]<>
						IntegerString[IntegerPart[100*dx],10,2]<>".txt",dx,memory];]]];
		]]]]]

calcAnh[fNames_List,masses_List,anhDir_String,nHess_Integer,nQ_Integer,dQ_,nModes_]:=
	Module[{mMat,mVec,hessData,hessCalc,n=Length[masses],s=Length[fNames],hessConv,
			mwh,vals,modes,freqs,phiTS,mp,tmp,mHessDat,pHessDat,mHessCalc,pHessCalc,
			pHessConv,mHessConv,Qsteps,phi3Calc,phi4Calc,phi3Conv,phi4Conv,i,j,k,l,
			cutOff=30/2.1947463*^5,fermi2,fermi3,anh,G0},
		hessData=Table[readInTxt[fNames[[i]],n,0][[2]],{i,s}];
		hessCalc=Table[If[i>=j,adDiff[hessData[[All,i,j]],2],
			Table[0,{s},{s}]],{i,3*n},{j,3*n}];
		hessConv=hessCalc[[All,All,s,s]];
		hessConv=hessConv+Transpose[hessConv]-DiagonalMatrix[Diagonal[hessConv]];
		mVec=Table[masses[[Floor[(i-1)/3]+1]],{i,3*n}];
		mMat=DiagonalMatrix[1/Sqrt[mVec]];
		mwh=mMat.hessConv.mMat;
		{vals,modes}=Eigensystem[mwh,nModes];
		mp=Convert[Dalton/ElectronMass,1];
		freqs=Sqrt[vals/mp];
		tmp=Position[Re[freqs],0.][[1]];
		freqs=Join[freqs[[tmp]],Delete[freqs,tmp]];
		modes=Join[modes[[tmp]],Delete[modes,tmp]];
		
		mHessDat=Table[readInTxt[anhDir<>"/"<>IntegerString[i]<>"m"<>
			IntegerString[{5,1,2,3,4,6,7,8,9,10,11,12}[[j]]]<>IntegerString[IntegerPart[100*0.01*2^(k-1)],10,2]<>
			".txt",n,0][[2]],{j,nModes},{i,nQ},{k,nHess}];
		pHessDat=Table[readInTxt[anhDir<>"/"<>IntegerString[i]<>"p"<>
			IntegerString[{5,1,2,3,4,6,7,8,9,10,11,12}[[j]]]<>IntegerString[IntegerPart[100*0.01*2^(k-1)],10,2]<>
			".txt",n,0][[2]],{j,nModes},{i,nQ},{k,nHess}];
		mHessCalc=Table[If[k>=l,adDiff[mHessDat[[j,i,All,k,l]],2],
			Table[0,{nHess},{nHess}]],{j,nModes},{i,nQ},{k,3*n},{l,3*n}];
		pHessCalc=Table[If[k>=l,adDiff[pHessDat[[j,i,All,k,l]],2],
			Table[0,{nHess},{nHess}]],{j,nModes},{i,nQ},{k,3*n},{l,3*n}];
		mHessConv=Table[If[k>=l,mHessCalc[[j,i,k,l,nHess,nHess]],mHessCalc[[j,i,l,k,nHess,nHess]]],
			{j,nModes},{i,nQ},{k,3*n},{l,3*n}];
		pHessConv=Table[If[k>=l,pHessCalc[[j,i,k,l,nHess,nHess]],pHessCalc[[j,i,l,k,nHess,nHess]]],
			{j,nModes},{i,nQ},{k,3*n},{l,3*n}];
		
		Qsteps=dQ/Sqrt[Abs[freqs]*mp];
		phiTS=modes.mMat.hessConv.mMat.Transpose[modes];

		phi3Calc=Table[0,{nModes},{nModes},{nModes},{3},{nQ},{nQ}];
		For[i=1,i<=nModes,i++,
		For[j=1,j<=nModes,j++,
		For[k=1,k<=nModes,k++,
			phi3Calc[[i,j,k,1]]=adDiff[Table[((modes.mMat.pHessConv[[i,l]].mMat.Transpose[modes])[[j,k]]-
				(modes.mMat.mHessConv[[i,l]].mMat.Transpose[modes])[[j,k]])/(2*Qsteps[[i]]*2^(l-1))/mp^(3/2),
				{l,nQ}],2];
			phi3Calc[[i,j,k,2]]=adDiff[Table[((modes.mMat.pHessConv[[j,l]].mMat.Transpose[modes])[[k,i]]-
				(modes.mMat.mHessConv[[j,l]].mMat.Transpose[modes])[[k,i]])/(2*Qsteps[[j]]*2^(l-1))/mp^(3/2),
				{l,nQ}],2];
			phi3Calc[[i,j,k,3]]=adDiff[Table[((modes.mMat.pHessConv[[k,l]].mMat.Transpose[modes])[[i,j]]-
				(modes.mMat.mHessConv[[k,l]].mMat.Transpose[modes])[[i,j]])/(2*Qsteps[[k]]*2^(l-1))/mp^(3/2),
				{l,nQ}],2];]]];
		phi3Conv=Table[Mean[phi3Calc[[i,j,k,{1,2,3},nQ,nQ]]],{i,nModes},{j,nModes},{k,nModes}];

		phi4Calc=Table[0,{nModes},{nModes},{2},{nQ},{nQ}];
		For[i=1,i<=nModes,i++,
		For[k=1,k<=nModes,k++,
			phi4Calc[[i,k,1]]=adDiff[Table[((modes.mMat.pHessConv[[k,l]].mMat.Transpose[modes])[[i,i]]+
				(modes.mMat.mHessConv[[k,l]].mMat.Transpose[modes])[[i,i]]-2*phiTS[[i,i]])/(Qsteps[[k]]*
				2^(l-1))^2/mp^2,{l,nQ}],2];
			phi4Calc[[i,k,2]]=adDiff[Table[((modes.mMat.pHessConv[[i,l]].mMat.Transpose[modes])[[k,k]]+
				(modes.mMat.mHessConv[[i,l]].mMat.Transpose[modes])[[k,k]]-2*phiTS[[k,k]])/(Qsteps[[i]]*
				2^(l-1))^2/mp^2,{l,nQ}],2];]];
		phi4Conv=Table[Mean[phi4Calc[[i,k,{1,2},nQ,nQ]]],{i,nModes},{k,nModes}];
		
		fermi2=Table[And[Im[freqs[[k]]]==0,Im[freqs[[l]]]==0,Abs[(phi3Conv[[k,k,l]]/freqs[[k]]/Sqrt[freqs[[l]]])^4/
			(2*freqs[[k]]-freqs[[l]])^3/256]>cutOff],{k,nModes},{l,nModes}];
		fermi3=Table[And[Im[freqs[[k]]]==0,Im[freqs[[l]]]==0,Im[freqs[[m]]]==0,Abs[phi3Conv[[k,l,m]]^4/freqs[[k]]^2/
			freqs[[l]]^2/freqs[[m]]^2/(freqs[[k]]-freqs[[l]]+freqs[[m]])^3/64]>cutOff],{k,nModes},{l,nModes},{m,nModes}];
		
		anh=Table[0,{nModes},{nModes}];
		For[k=1,k<=nModes,k++,
			anh[[k,k]]=1/16/freqs[[k]]^2*(phi4Conv[[k,k]]-Sum[phi3Conv[[k,k,l]]^2/2/freqs[[l]]*(1/(2*freqs[[k]]+freqs[[l]])-
				If[fermi2[[k,l]],0,1/(2*freqs[[k]]-freqs[[l]])]+4/freqs[[l]]),{l,nModes}]);];
		For[k=1,k<=nModes,k++,
		For[l=1,l<=nModes,l++,
			If[k!=l,anh[[k,l]]=1/4/freqs[[k]]/freqs[[l]]*(phi4Conv[[k,l]]-Sum[phi3Conv[[k,k,m]]*phi3Conv[[l,l,m]]/freqs[[m]]^2,
				{m,nModes}]+Sum[phi3Conv[[k,l,m]]^2/2/freqs[[m]]*(-1/(freqs[[k]]+freqs[[l]]+freqs[[m]])+If[fermi3[[k,m,l]],
				0,1/(freqs[[k]]+freqs[[l]]-freqs[[m]])]-If[fermi3[[k,l,m]],0,1/(freqs[[k]]-freqs[[l]]+freqs[[m]])]+
				If[fermi3[[l,k,m]],0,1/(freqs[[k]]-freqs[[l]]-freqs[[m]])]),{m,nModes}])];]];
		
		G0=Sum[phi4Conv[[i,i]]/freqs[[i]]^2,{i,nModes}]/64-7/576*Sum[phi3Conv[[i,i,i]]^2/freqs[[i]]^4,{i,nModes}]+
			3/64*Sum[If[k==l,0,phi3Conv[[k,l,l]]^2/freqs[[l]]^2*If[fermi2[[l,k]],0,1/(2*freqs[[l]]-freqs[[k]])]*
			1/(2*freqs[[l]]+freqs[[k]])],{k,nModes},{l,nModes}]-Sum[If[k<l<m&&!(fermi3[[k,m,l]]||fermi3[[k,l,m]]),
			phi3Conv[[k,l,m]]^2/((freqs[[k]]+freqs[[l]])^2-freqs[[m]]^2)/((freqs[[k]]-freqs[[l]])^2-freqs[[m]]^2),0],
			{m,nModes},{l,nModes},{k,nModes}]/4;
		{anh,G0}
		]

anh1D[en_List,dQ_,freq_]:=
	Module[{thirdD,fourthD,xFF,G0,nE=Length[en],deltaQ,tsI,nAD},
		deltaQ=dQ/Sqrt[Abs[freq]];
		tsI=Floor[nE/2]+1;
		nAD=nE-tsI-1;
		thirdD=adDiff[Table[(en[[tsI+i+1]]-2*en[[tsI+i]]+2*en[[tsI-i]]-en[[tsI-i-1]])/
			2/(deltaQ*2^(i-1))^3,{i,nAD}],2][[nAD,nAD]];
		fourthD=adDiff[Table[(en[[tsI+i+1]]-4*en[[tsI+i]]+6*en[[tsI]]-4*en[[tsI-i]]+
			en[[tsI-i-1]])/(deltaQ*2^(i-1))^4,{i,nAD}],2][[nAD,nAD]];
		xFF=1/16/freq^2*(fourthD-5*thirdD^2/3/freq^2);
		G0=fourthD/freq^2/64-5/576*thirdD^2/freq^4;
		{xFF,G0}]/;(Length[en]/2-Floor[Length[en]/2])==0.5

End[]
EndPackage[]
