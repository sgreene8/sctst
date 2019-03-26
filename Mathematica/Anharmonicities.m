(* ::Package:: *)

BeginPackage["Anharmonicities`"]
calc2DAnh::usage="calc2DAnh[deltaQ,masses,projVec] uses the MOLPRO output files ts.txt, p1.txt, m1.txt, p2.txt, and m2.txt in the current directory to calculate a matrix of anharmonicity constants using rectilinear projection. The step size, in Sqrt[Da]a0, for each mode is inputted in the deltaQ vector. Alternatively, 'CFOUR' can be inputted for deltaQ to use the default CFOUR step sizes. nAtom is the total number of atoms in the system. masses is a list of masses corresponding to each atom in the system. projVec defines the 2 bonds being treated in this RD method, as in RectilProj.m."
curvilRDAnh::usage="curvilRDAnh[deltaQ,internalcoord,masses,projVec,mNames,pNames] calculates a matrix of anharmonicity constants using curvilinear projection, using the ts.txt file and files with names specified in pNames and mNames in the current directory. The step size, in Sqrt[Da]a0, for each mode is inputted in the deltaQ vector. Alternatively, 'CFOUR' can be inputted for deltaQ to use the default CFOUR step sizes. nAtom is the total number of atoms in the system. masses is a list of masses corresponding to each atom in the system. projVec defines the 2 curvilinear coordinates being treated in this RD method, as in CurvilProj.m."
calcFDAnh::usage="calcFDAnh[deltaQ,masses,mNames,pNames] calculates a matrix of anharmonicity constants for all 3n-6 vibrational degrees of freedom, using the ts.txt file and files with names specified in pNames and mNames in the current directory. The step size, in Sqrt[Da]a0, for each mode is inputted in the deltaQ vector. Alternatively, 'CFOUR' can be inputted for deltaQ to use the default CFOUR step sizes. nAtom is the total number of atoms in the system. masses is a list of masses corresponding to each atom in the system."


Begin["`Private`"]
Needs["readInTxt`"]
Needs["RectilProj`"]
Needs["CurvilProj`"]

calc2DAnh[dQ_,mVec_List,pV_]:=
	Module[{nAtom=Length[mVec],tsHess,tsGeom,mMat,p,pmat,L,phiTS,
			mp=1.660538922*^-27/9.1093897*^-31,tsFreqAU,phiTen,hess,fNames,thirdD,fourthD,anh,
			deltaQ,G0},
		{tsGeom,tsHess}=readInTxt["ts.txt",nAtom,1][[{1,2}]];
		mMat=Table[0,{3 nAtom},{3 nAtom}];
		Do[mMat[[i,i]]=mVec[[Floor[(i-1)/3]+1]]^-0.5,{i,1,3 nAtom}];
		p=Table[0,{2}];
		Do[p[[i]]=RectilBond[tsGeom,mVec,pV[[i]]],{i,1,2}];
        pmat=Transpose[Orthogonalize[p]].Orthogonalize[p];
		L=Eigenvectors[pmat.mMat.tsHess.mMat.pmat,2];
		(*This matrix has units of a.u.*)
		phiTS=L.mMat.tsHess.mMat.Transpose[L]/mp;
		tsFreqAU=Sqrt[Diagonal[phiTS]];
		
		(*Phi tensor is {{phip1, phim1}, {phip2, phim2}}*)
		phiTen=Table[0,{2},{2},{2},{2}];
		fNames={{"p1.txt","m1.txt"},{"p2.txt","m2.txt"}};
		Do[hess=readInTxt[fNames[[i,j]],nAtom,1][[2]];
			phiTen[[i,j]]=L.mMat.hess.mMat.Transpose[L]/mp,{i,2},{j,2}];
		deltaQ=If[And[StringQ[dQ],dQ=="CFOUR"],0.05/Sqrt[Abs[tsFreqAU]mp],dQ];
		thirdD=Table[Mean[{(phiTen[[i,1,j,k]]-phiTen[[i,2,j,k]])/(2deltaQ[[i]]Sqrt[mp]),
			(phiTen[[j,1,k,i]]-phiTen[[j,2,k,i]])/(2deltaQ[[j]]Sqrt[mp]),(phiTen[[k,1,i,j]]-
			phiTen[[k,2,i,j]])/(2deltaQ[[k]]Sqrt[mp])}],{i,2},{j,2},{k,2}];
		fourthD=Table[If[i==j,Mean[{(phiTen[[k,1,i,j]]+phiTen[[k,2,i,j]]-2phiTS[[i,j]])/mp/deltaQ[[k]]^2,
			(phiTen[[i,1,k,k]]+phiTen[[i,2,k,k]]-2phiTS[[k,k]])/deltaQ[[i]]^2/mp}],
			(phiTen[[k,1,i,j]]+phiTen[[k,2,i,j]]-2phiTS[[i,j]])/deltaQ[[k]]^2/mp],{i,2},{j,2},{k,2}];
		anh=Table[0,{2},{2}];
		Do[anh[[k,k]]=1/(16*tsFreqAU[[k]]^2)*(fourthD[[k,k,k]]-Sum[thirdD[[k,k,i]]^2*
			(8tsFreqAU[[k]]^2-3tsFreqAU[[i]]^2)/(tsFreqAU[[i]]^2*(4tsFreqAU[[k]]^2-
			tsFreqAU[[i]]^2)),{i,2}]),{k,2}];
		anh[[1,2]]=1/(4tsFreqAU[[1]]tsFreqAU[[2]])*(fourthD[[1,1,2]]-Sum[thirdD[[1,1,m]]
			thirdD[[2,2,m]]/tsFreqAU[[m]]^2,{m,2}]+Sum[2thirdD[[1,2,m]]^2*(tsFreqAU[[1]]^2+
			tsFreqAU[[2]]^2-tsFreqAU[[m]]^2)/(((tsFreqAU[[1]]+tsFreqAU[[2]])^2-tsFreqAU[[m]]^2)*
			((tsFreqAU[[1]]-tsFreqAU[[2]])^2-tsFreqAU[[m]]^2)),{m,2}]);
		anh[[2,1]]=anh[[1,2]];
		
		G0=Sum[fourthD[[i,i,i]]/tsFreqAU[[i]]^2,{i,2}]/64-Sum[thirdD[[i,i,i]]^2/tsFreqAU[[i]]^4,{i,2}]*7/576+
			Sum[If[k==l,0,thirdD[[k,l,l]]^2/(4tsFreqAU[[l]]^2-tsFreqAU[[k]]^2)/tsFreqAU[[l]]^2],{k,2},{l,2}]*3/64
			-Sum[thirdD[[k,l,m]]^2/((tsFreqAU[[k]]+tsFreqAU[[l]])^2-tsFreqAU[[m]]^2)/
			((tsFreqAU[[k]]-tsFreqAU[[l]])^2-tsFreqAU[[m]]^2),{m,2},{l,1,m},{k,1,l}]/4;
		
		{anh,Re[G0]}]

curvilRDAnh[dQ_,int_List,m_List,pV_List,mNames_,pNames_]:=
	Module[{nAtom=Length[m],coord,hess,grid,l,li=Length[int],lp=Length[pV],mV,mMat,freq,cartvec,phits,phim,phip,
phi3,phi4,mp=1.660538922*10^(-27)/(9.1093897*10^(-31)),anh,dd,G0},
{coord,hess,grid}=readInTxt["ts.txt",nAtom,0][[{1,2,3}]];
l=3Length[coord];
mV=Table[0,{l}];
Do[mV[[i]]=m[[Floor[(i-1)/3]+1]],{i,1,l}];
mMat=DiagonalMatrix[1/Sqrt[mV]];
cartvec=curvilProjRD[coord,int,hess,grid,m,pV,{0,0},"vec"];
phits=cartvec.(mMat.hess.mMat).Transpose[cartvec];

freq=Sqrt[Diagonal[phits]/mp];
dd=If[And[StringQ[dQ],dQ=="CFOUR"],0.05/Sqrt[Abs[freq]mp],dQ];
Print[dd];
phim=phip=Table[0,{lp}];
Do[phim[[i]]=cartvec.(mMat.readInTxt[mNames[[i]],nAtom,0][[2]].mMat).Transpose[cartvec],{i,1,lp}];
Do[phip[[i]]=cartvec.(mMat.readInTxt[pNames[[i]],nAtom,0][[2]].mMat).Transpose[cartvec],{i,1,lp}];
phi3=Table[0,{lp},{lp},{lp}];
Do[phi3[[i,j,k]]=1/6((phip[[i]][[j,k]]-phim[[i]][[j,k]])/dd[[i]]+(phip[[j]][[i,k]]-phim[[j]][[i,k]])/dd[[j]]+(phip[[k]][[i,j]]-phim[[k]][[i,j]])/dd[[k]]),{i,1,lp},{j,1,lp},{k,1,lp}];
phi3=phi3/Sqrt[mp]^3;
phi4=Table[0,{lp},{lp},{lp},{lp}];

Do[phi4[[i,i,k,k]]=1/2((phip[[k]][[i,i]]+phim[[k]][[i,i]]-2phits[[i,i]])/dd[[k]]^2+(phip[[i]][[k,k]]+phim[[i]][[k,k]]-2phits[[k,k]])/dd[[i]]^2),{i,1,lp},{k,1,lp}];
phi4=phi4/mp^2;
anh=Table[0,{lp},{lp}];
Do[anh[[k,i]]=If[k<i,0,1/(4freq[[k]]freq[[i]])(phi4[[k,k,i,i]]-Sum[phi3[[k,k,j]]phi3[[i,i,j]]/freq[[j]]^2,{j,1,lp}]+Sum[2phi3[[k,i,j]]^2(freq[[k]]^2+freq[[i]]^2-freq[[j]]^2)/(((freq[[k]]+freq[[i]])^2-freq[[j]]^2)((freq[[k]]-freq[[i]])^2-freq[[j]]^2)),{j,1,lp}])],{k,1,lp},{i,1,lp}];
Do[anh[[k,k]]=1/(16freq[[k]]^2)(phi4[[k,k,k,k]]-Sum[phi3[[k,k,i]]^2(8freq[[k]]^2-3freq[[i]]^2)/(freq[[i]]^2(4freq[[k]]^2-freq[[i]]^2)),{i,1,lp}]),{k,1,lp}];
anh=Transpose[anh]+anh-DiagonalMatrix[Diagonal[anh]];

G0=Sum[phi4[[k,k,k,k]]/freq[[k]]^2,{k,lp}]/64-Sum[phi3[[k,k,k]]^2/freq[[k]]^4,{k,lp}]*7/576+
Sum[If[k==i,0,phi3[[k,i,i]]^2/(4freq[[i]]^2-freq[[k]]^2)/freq[[i]]^2],{k,lp},{i,lp}]*3/64-
Sum[phi3[[k,i,a]]^2/((freq[[k]]+freq[[i]])^2-freq[[a]]^2)/((freq[[k]]-freq[[i]])^2-freq[[a]]^2),{a,lp},{i,1,a},{k,1,i}]/4;

{anh,Re[G0]}]/;Length[mNames]==Length[pNames]==Length[pV]

calcFDAnh[dQ_,m_List,mNames_List,pNames_List]:=
	Module[{nAtom=Length[m],l,mV,mMat,coord,hess,freq,cartvec,phits,phim,phip,phi3,phi4,
mp=1.660538922*10^(-27)/(9.1093897*10^(-31)),anh,dd},
	l=3nAtom-6;
	mV=Table[0,{l}];
	Do[mV[[i]]=m[[Floor[(i-1)/3]+1]],{i,1,l+6}];
	mMat=DiagonalMatrix[1/Sqrt[mV]];
	{coord,hess}=readInTxt["ts.txt",nAtom,0][[{1,2}]];
	cartvec=RectilProj[coord,hess,m,{{1,2},{2,3}},"pre"][[2]];
	cartvec=Table[cartvec[[i]],{i,l}];
	phits=cartvec.(mMat.hess.mMat).Transpose[cartvec];
	freq=Sqrt[Diagonal[phits]/mp];
	dd=If[And[StringQ[dQ],dQ=="CFOUR"],0.05/Sqrt[Abs[freq]mp],dQ];
	phim=phip=Table[0,{l-6}];
	Do[phim[[i]]=cartvec.(mMat.readInTxt[mNames[[i]],nAtom,0][[2]].mMat).Transpose[cartvec],{i,1,l}];
	Do[phip[[i]]=cartvec.(mMat.readInTxt[pNames[[i]],nAtom,0][[2]].mMat).Transpose[cartvec],{i,1,l}];
	phi3=Table[0,{l},{l},{l}];
	Do[phi3[[i,j,k]]=1/6((phip[[i]][[j,k]]-phim[[i]][[j,k]])/dd[[i]]+(phip[[j]][[i,k]]-phim[[j]][[i,k]])/dd[[j]]+(phip[[k]][[i,j]]-phim[[k]][[i,j]])/dd[[k]]),{i,1,l},{j,1,l},{k,1,l}];
	phi3=phi3/Sqrt[mp]^3;
	phi4=Table[0,{l},{l},{l},{l}];
	Do[phi4[[i,i,k,k]]=1/2((phip[[k]][[i,i]]+phim[[k]][[i,i]]-2phits[[i,i]])/dd[[k]]^2+(phip[[i]][[k,k]]+phim[[i]][[k,k]]-2phits[[k,k]])/dd[[i]]^2),{i,1,l},{k,1,l}];
	Do[anh[[k,i]]=If[k<i,0,1/(4freq[[k]]freq[[i]])(phi4[[k,k,i,i]]-Sum[phi3[[k,k,j]]phi3[[i,i,j]]/freq[[j]]^2,{j,1,l}]+Sum[2phi3[[k,i,j]]^2(freq[[k]]^2+freq[[i]]^2-freq[[j]]^2)/(((freq[[k]]+freq[[i]])^2-freq[[j]]^2)((freq[[k]]-freq[[i]])^2-freq[[j]]^2)),{j,1,l}])],{k,1,l},{i,1,l}];
	Do[anh[[k,k]]=1/(16freq[[k]]^2)(phi4[[k,k,k,k]]-Sum[phi3[[k,k,i]]^2(8freq[[k]]^2-3freq[[i]]^2)/(freq[[i]]^2(4freq[[k]]^2-freq[[i]]^2)),{i,1,l}]),{k,1,l}];
	anh=Transpose[anh]+anh-DiagonalMatrix[Diagonal[anh]];
	anh]/;Length[mNames]==Length[pNames]==(3Length[m]-6)
End[]
EndPackage[]
