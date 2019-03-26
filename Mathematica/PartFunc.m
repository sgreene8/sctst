(* ::Package:: *)
(* Contributed by Xiao Shan *)
(* ::Input:: *)
(**)


BeginPackage["PartFunc`"]
PartFunc::usage="This package contains five new keywords to produce the non electronic partition function as a function of temperature for molecules and atoms. The new keywords are rotN, nonElecPartMol, nonElecPartAtom, nonElecPartTSRecP and nonElecPartTSCurP. It includes both rectilinear and curvilinear projections to work out the spectator mode partition function for a transition state. The RectilProj.m and CurvilProj.m packages must have been installed."
rotN::usage="rotN[string] produces rotational number of a certain point group, detailed discussion can be found in Theor. Chem. Account 118 (2007) 813"
nonElecPartMol::usage="nonElecPartMol[coord,hess,m,pg,t] produces the partition function, including rotational, translational and vibrational ones, of a molecule (note: transition state before projection is included) w.r.t. temperature, t. The coordinate, hessian matrix, mass vector, point group are required as input."
nonElecPartMolIR::usage="nonElecPartMol[coord,hess,m,pg,t,rotAtoms,rotAxis] produces the partition function, including rotational, translational, vibrational, and internal rotation, of a molecule (note: transition state before projection is included) w.r.t. temperature, t. The coordinate, hessian matrix, mass vector, point group are required as input. The molecule is assumed to have 1 internal rotation of symmetry 3 with the atoms specified in rotAtoms[[1]] rotating relative to those in rotAtoms[[2]], with the rotation axis defined by the 2 atoms in rotAxis."
nonElecPartAtom::usage="nonElecPartAtom[m,t] produces the partition function of atom w.r.t. temperature, t. The atomic mass in a.u. is required as input."
nonElecPartTSRP::usage="nonElecPartTSRP[coord,hess,m,projVec,pg,t] produces the partition function, including rotational, translational and spectator mode (according to rectilinear projection method) vibrational ones, of a transition state structure w.r.t. to temperature, t. The coordinate, hessian matrix, mass vector, projection vector and point group are required as input."
nonElecPartTSRPIR::usage="nonElecPartTSRPIR[coord,hess,m,projVec,pg,t,rotAtoms,rotAxis] produces the partition function, including rotational, translational and spectator mode (according to rectilinear projection method) vibrational, and internal rotation, of a transition state structure w.r.t. to temperature, t. The coordinate, hessian matrix, mass vector, projection vector and point group are required as input. The molecule is assumed to have 1 internal rotation of symmetry 3 with the atoms specified in rotAtoms[[1]] rotating relative to those in rotAtoms[[2]], with the rotation axis defined by the 2 atoms in rotAxis."
nonElecPartTSCP::usage="nonElecPartTSCP[coord,hess,grad,m,internalCoordVec,projVec,pg,t] produces the partition function, including rotational, translational and spectator mode (according to curvil projection method) vibrational ones, of a transition state structure w.r.t. to temperature, t. The coordinate, hessian matrix, gradient, mass vector, internal coordinate vector, projection vector and point group are required as input."
nonElecPartTSCPIR::usage="nonElecPartTSCPIR[coord,hess,grad,m,internalCoordVec,projVec,pg,t,rotAtoms,rotAxis] produces the partition function, including rotational, translational and spectator mode (according to curvil projection method) vibrational, and internal rotation, of a transition state structure w.r.t. to temperature, t. The coordinate, hessian matrix, gradient, mass vector, internal coordinate vector, projection vector and point group are required as input. The molecule is assumed to have 1 internal rotation of symmetry 3 with the atoms specified in rotAtoms[[1]] rotating relative to those in rotAtoms[[2]], with the rotation axis defined by the 2 atoms in rotAxis."
nonElecnonVibPartTS::usage="nonElecnonVibPartTS[coord,m,pg,t] produces only the product of the translational and rotational partition functions for the transition state. The coordinate, mass vector, and point group are required as input."
pointLineDist::usage="pointLineDist[l1,l2,p] calculates the distance between a point p and the line that passes through points l1 and l2 according to the formula from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html."
qhr::usage="qhr[coord,m,rotAtoms,rotAxis,freq,t] calculates the partition fxn for a hindered rotor. rotAtoms is a 2-D list of the 2 groups of atoms that are rotating (their indices in the inputted coordinates), and rotAxis defines the bond about which they rotate. freq is the harmonic frequency associated with the rotation, and t is the temperature. Assumes 3-fold symmetry of rotation."
qvib::usage="qhr[freqs,t] calculates a vibrational partition function for a given set of frequencies."
qrot::usage="qrot[I,pg,t] calculates a rotational partition function given a vector of moments of inertia (Da a0^2). If the vector only 2 nonzero elements, the molecule is assumed to be linear."

Begin["`Private`"]
Needs["CurvilProj`"];
Needs["RectilProj`"];

rotN[x_String]:=
 Module[{n,t,i},
        t=Transpose[{{"c1","cs","c2","c2v","c3v","c00v","d2h","d3h","d5h","d00h","d3d","td","oh"},
                     {1,1,2,2,3,1,4,6,10,2,6,12,24}}];
        n=0;
        For[i=1,i<=Length[t],i++,If[x==t[[i,1]],n=t[[i,2]],n=n]];
       n]

pointLineDist[l1_List,l2_List,p_List]:=
	Module[{},Norm[Cross[p-l1,p-l2]]/Norm[l1-l2]]/;Length[l1]==Length[l2]==Length[p]==3


nonElecPartMol[coord_,hess_,m_List,pg_String,t_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),mV,l=Length[coord],freq,qtrans,intmat,int1},
        mV=Table[0,{3l}];
        Do[mV[[i]]=m[[Floor[(i-1)/3]+1]],{i,1,3l}];
        freq=Sqrt[Eigenvalues[DiagonalMatrix[1/Sqrt[mV]].hess.DiagonalMatrix[1/Sqrt[mV]]]]*5.14048695*10^3;
        Do[freq[[i]]=If[Abs[freq[[i]]]<=10,0,If[Re[freq[[i]]]==0,0,freq[[i]]]],{i,1,3l}];
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
        Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,k]]^2,{k,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int1=Eigenvalues[intmat];
        qvib[freq,t]qtrans[t]qrot[int1,pg,t]
        ]/;3Length[coord]==3Length[m]==Dimensions[hess][[1]]==Dimensions[hess][[2]]

qvib[freqs_List,t_]:=
	Module[{hbar=1.054571*^-34,kb=1.380650*^-23,ccm=2.997925*10^10,l=Length[freqs]},
		Product[If[Re[freqs[[i]]]==0,1,1/(1-Exp[-2Pi hbar*ccm/(kb t)*freqs[[i]]])],{i,l}]
]

qrot[moments_List,pg_String,t_]:=
	Module[{autom=5.291772*^-11,hbar=1.054571*^-34,kb=1.380650*^-23,autokg=1.660539*^-27,rotn},
		(*The if statement accounts for linear molecules.*)
		Sqrt[Pi*t^3/Product[hbar^2/(2*moments[[n]] autom^2 autokg kb),{n,2}]/
			If[moments[[3]]==0,Pi*t,hbar^2/(2*moments[[3]] autom^2 autokg kb)]]/rotN[pg]
]

qhr[coord_,m_List,rotAtoms_,rotAxis_,freq_,t_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),autom=5.291772*10^(-11),
		i1,i2,ieff,qi,qho,qfr},
		(*See Chuang & Truhlar (2000)*)
		qho[x_]:=Exp[-2Pi hbar*ccm/(kb x)*freq/2]/(1-Exp[-2Pi hbar*ccm/(kb x)*freq]);
		(*qho[x_]:=1/(1-Exp[-2Pi hbar*ccm/(kb x)*freq]);*)
		i1=Sum[m[[rotAtoms[[1,i]]]]pointLineDist[coord[[rotAxis[[1]]]],coord[[rotAxis[[2]]]],coord[[rotAtoms[[1,i]]]]]^2,{i,Length[rotAtoms[[1]]]}];
		i2=Sum[m[[rotAtoms[[2,i]]]]pointLineDist[coord[[rotAxis[[1]]]],coord[[rotAxis[[2]]]],coord[[rotAtoms[[2,i]]]]]^2,{i,Length[rotAtoms[[2]]]}];
		ieff=i1*i2/(i1+i2);
		qfr[x_]:=(2Pi*ieff*autokg*autom^2*kb*x)^0.5/(hbar 3);
		qi[x_]:=kb*x/(2Pi*hbar*ccm*freq);
		qho[t]*Tanh[qfr[t]/qi[t]]
]/;Length[coord]==Length[m]

nonElecPartMolIR[coord_,hess_,m_List,pg_String,t_,rotAtoms_,rotAxis_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),
		mV,l=Length[coord],freq,qtrans,intmat,int1,i1,i2,ieff,qfr,qho,qi},
        mV=Table[0,{3l}];
        Do[mV[[i]]=m[[Floor[(i-1)/3]+1]],{i,1,3l}];
        freq=Sqrt[Eigenvalues[DiagonalMatrix[1/Sqrt[mV]].hess.DiagonalMatrix[1/Sqrt[mV]],3l-6]]*5.14048695*10^3;
		freq=Re[freq];
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
        Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,k]]^2,{k,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int1=Eigenvalues[intmat];
        qvib[Delete[freq,-1],t]qtrans[t]qrot[int1,pg,t]*qhr[coord,m,rotAtoms,rotAxis,freq[[3l-6]],t]
        ]/;And[3Length[coord]==3Length[m]==Dimensions[hess][[1]]==Dimensions[hess][[2]],Length[rotAtoms]==Length[rotAxis]==2]

nonElecPartAtom[m_,t_]:=
 Module[{autokg=1.660539*10^(-27),kb=1.380650*10^(-23),h=6.626*10^(-34)},
		Sqrt[2Pi m autokg kb t/h^2]^3]

nonElecPartTSRP[coord_,hess_,m_List,pV_,pg_String,t_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),l=Length[coord],freq,qtrans,intmat,int1},
        freq=RectilProj[coord,hess,m,pV,"post"][[1]];
		Do[freq[[i]]=If[Abs[freq[[i]]]<=10,0,If[Re[freq[[i]]]==0,0,freq[[i]]]],{i,1,3l}];
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
		Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,q]]^2,{q,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int1=Eigenvalues[intmat];
        qvib[freq,t]qtrans[t]qrot[int1,pg,t]
        ]/;3Length[coord]==3Length[m]==Dimensions[hess][[1]]==Dimensions[hess][[2]]

nonElecPartTSRPIR[coord_,hess_,m_List,pV_,pg_String,t_,rotAtoms_,rotAxis_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),
		l=Length[coord],freq,qtrans,intmat,int1,i1,i2,ieff,qfr,qho,qi,qhr},
        freq=RectilProj[coord,hess,m,pV,"post"][[1]];
		Do[freq[[i]]=If[Abs[freq[[i]]]<=10,0,If[Re[freq[[i]]]==0,0,freq[[i]]]],{i,1,3l}];
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
		Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,q]]^2,{q,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int1=Eigenvalues[intmat];
        qvib[freq[[1;;3l-9]],t]qtrans[t]qrot[int1,pg,t]*qhr[coord,m,rotAtoms,rotAxis,freq[[3l-8]],t]
        ]/;3Length[coord]==3Length[m]==Dimensions[hess][[1]]==Dimensions[hess][[2]]&&Length[rotAtoms]==Length[rotAxis]==2

nonElecPartTSCP[coord_,hess_,grad_,m_List,intcoord_,pV_,pg_String,t_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),l=Length[coord],freq,qtrans,intmat,int},
        freq=curvilProj[coord,intcoord,hess,grad,m,pV,"post"][[2]];
        Do[freq[[i]]=If[Abs[freq[[i]]]<=10,0,If[Re[freq[[i]]]==0,0,freq[[i]]]],{i,1,Length[freq]}];
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
        Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,k]]^2,{k,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int=Eigenvalues[intmat];
        qvib[freq,t]qtrans[t]qrot[int,pg,t]
        ]/;3Length[coord]==3Length[m]==Dimensions[hess][[1]]==Dimensions[hess][[2]]==Length[grad]

nonElecPartTSCPIR[coord_,hess_,grad_,m_List,intcoord_,pV_,pg_String,t_,rotAtoms_,rotAxis_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),
		l=Length[coord],freq,qtrans,intmat,int,i1,i2,ieff,qfr,qho,qi,qhr},
        freq=curvilProj[coord,intcoord,hess,grad,m,pV,"post"][[2]];
        Do[freq[[i]]=If[Abs[freq[[i]]]<=10,0,If[Re[freq[[i]]]==0,0,freq[[i]]]],{i,1,Length[freq]}];
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
        Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,k]]^2,{k,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int=Eigenvalues[intmat];
        qvib[freq[[1;;3l-9]],t]qtrans[t]qrot[int,pg,t]*qhr[coord,m,rotAtoms,rotAxis,freq[[3l-8]],t]
        ]/;And[3Length[coord]==3Length[m]==Dimensions[hess][[1]]==Dimensions[hess][[2]]==Length[grad],Length[rotAtoms]==Length[rotAxis]==2]

nonElecnonVibPartTS[coord_,m_List,pg_String,t_]:=
 Module[{hbar=1.054571*10^(-34),kb=1.380650*10^(-23),ccm=2.997925*10^10,autokg=1.660539*10^(-27),l=Length[coord],qtrans,intmat,int},
        qtrans[x_]:=(Sqrt[2Pi Sum[m[[i]],{i,1,l}]autokg kb x/(2Pi hbar)^2])^3;
        intmat=Table[0,{3},{3}];
        Do[intmat[[i,j]]=If[i==j,Sum[m[[n]](Sum[coord[[n,k]]^2,{k,1,3}]-coord[[n,i]]^2),{n,1,l}],-Sum[m[[n]]coord[[n,i]]coord[[n,j]],{n,1,l}]],{i,1,3},{j,1,3}];
        int=Eigenvalues[intmat];
Print[int];
        qtrans[t]qrot[int,pg,t]
        ]/;3Length[coord]==3Length[m]

End[]
EndPackage[]

