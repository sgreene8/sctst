(* ::Package:: *)
(* Contributed by Xiao Shan *)

BeginPackage["readInTxt`"]
readInTxt::usage="readInTxt[filename,n,ne] reads in a txt file, which consists the data from a molpro frequency calculation output, n is the total number of atoms in the system, ne is the total number of single point energies associated with the geometry. A typical output txt file must consists atomic coordinate, force constants, numerical gradients and energy."
readInTxt::usage="readInTxt[filename,n,ne] produces a table of 4 elements corresponding to the atomic coordinates, hessian matrix, numerical gradients and energy."
readInTxtAtom::usage="readInTxtAtom[filename,ne] produces a list of ne elements corresponding to the energies of an atom."

 
Begin["`Private`"]
readInTxt[filename_String,n_Integer,ne_Integer]:=Module[{file,t=Table[0,{4}],a},
 file=OpenRead[filename];
 Do[t[[1]]=Table[0,{n},{3}]];
 Do[t[[2]]=Table[0,{3n},{3n}]];
 Do[t[[3]]=Table[0,{3n}]];
 Do[t[[4]]=Table[0,{ne}]];
 Find[file,"Atomic Coordinates"];
 For[i=1,i<=n,i++,For[j=1,j<=3,j++,t[[1]][[i,j]]=Read[file,Number]]];
 Find[file,"Force Constants"];
 For[k=1,k<=Floor[3n/5]+1,k++,For[i=(k-1)*5+1,i<=3n,i++,For[j=(k-1)*5+1,j<=Min[i,Min[k*5,3n]],j++,
     t[[2]][[i,j]]=Read[file,Number]]]];
 Do[t[[2]]=t[[2]]+Transpose[t[[2]]]-DiagonalMatrix[Diagonal[t[[2]]]]];
 a=Find[file,"Numerical gradient"];
 If[a===EndOfFile,Find[file,"Numerical Gradient"]];
 If[a===EndOfFile,
 For[i=1,i<=3n,i++,t[[3]][[i]]=0];
 SetStreamPosition[file,1],
 For[i=1,i<=3n,i++,t[[3]][[i]]=Read[file,Number]]];
 Find[file,"Energy"];
 For[i=1,i<=ne,i++,t[[4]][[i]]=Read[file,Number]];
 Close[filename];
 t]

readInTxtAtom[filename_String,ne_Integer]:=
 Module[{file,t=Table[0,{ne}]},
        file=OpenRead[filename];
        Find[file,"Energy"];
        For[i=1,i<=ne,i++,t[[i]]=Read[file,Number]];
        Close[filename];
        t]
 
 
End[]
EndPackage[]
