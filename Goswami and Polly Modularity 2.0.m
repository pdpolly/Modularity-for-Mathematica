(* ::Package:: *)

(*  This function prints the version number for the installed verison of this package.  *)

ModularityVersion[]:=Print["Modularity for Mathematica 2.0\n(c) P. David Polly and A. Goswami, 20 Dec 2010\nUpdated 7 April 2018."];
ModularityVersion[];


(*  Modularity Package, written by P.David Polly and Anjali Goswami 2010.

This package performs a simple analysis of modularity for geometric morphometric landmark data. 

Three functions for landmarks are included, two of which calculate correlation matrices for multidimensional 
data (congruence coefficient and canonical correlation respectively) and one which uses those correlation 
matrices to perform a simple cluster analysis of modularity. An additional function, Magwene[],performs a graphical 
modeling analysis as described by Magwene (2001).The algorithms in the Magwene[] function do not work for singular 
correlation matrices,such as those derived from Procrustes superimposed landmarks.For further information see 
Goswami and Polly (2010). Functions for performing Procrustes analysis, matrix correlation analysis, Mantel's test, and 
subsampling analysis are also included.

Goswami,A. and P.D.Polly. 2010. Methods for studying morphological integration,modularity and covariance 
evolution.Pp.213-243 in J.Alroy and G.Hunt,eds., Quantitative Methods in Paleobiology. Paleontological 
Society Short Course, October 30th, 2010. The Paleontological Society Papers,Volume 16.

Magwene,P. 2001. New tools for studying integration and modularity. Evolution, 55:1734-1745.
*)



(* This functions does a Procrustes superimposition of landmark coordinates followed by orthogonal projection into tangent space.  The algorithm is the one 
    presented by Rohlf & Slice, 1990, Syst. Zool. 39: 40-59  with tangent space projection from Rohlf, 1999, J. Class. 16: 197-223.

   Usage:    Procrustes[data, n, k] 
             where data are the landmark coodiantes to be aligned, n is the number of landmarks, and k is the number of dimensions of each landmark (2 or 3).
             Aligned coordinates are returned with each shape in a single row, with the columns being the x,y (,z) coordinates of the n landmarks.

	Updated 21 February 2010 to fix problem with 3D rotations.
	Updated 7 April 2018 to be consistent with current versions of Mathematica.  
*)
Procrustes[data_,nlandmarks_,ndims_]:=Module[{l,II,PP,SS,x,y,u,w,v,hh,ResidSS,NewResidSS},
l=Partition[Partition[Flatten[data],ndims],nlandmarks];
II=IdentityMatrix[nlandmarks];
PP=Table[N[1/nlandmarks],{nlandmarks},{nlandmarks}];
SS=Table[N[Sqrt[Tr[(II-PP).l[[x]].Transpose[l[[x]]].(II-PP)]]],{x,Length[l]}];
l=Table[((II-PP).l[[x]])/SS[[x]],{x,Length[l]}];
y=Mean[l];
ResidSS=Plus@@Flatten[Table[(Flatten[l[[x]]]-Flatten[y])^2,{x,Length[l]}]];
While[True,
	For[x=1,x<=Length[l],x++,
		{u,w,v}=SingularValueDecomposition[Transpose[l[[x]]].y];
		hh=u.(w*Inverse[Abs[w]]).Transpose[v];
		l[[x]]=l[[x]].hh;

		];
	y=Mean[l];
    NewResidSS=Plus@@Flatten[Table[(Flatten[l[[x]]]-Flatten[y])^2,{x,Length[l]}]];
	If[Abs[NewResidSS-ResidSS]<0.0001,Break[]];
	ResidSS=NewResidSS;
];
l=Partition[PrincipalComponents[Partition[Flatten[l],ndims]],nlandmarks];
l=Partition[Flatten[l],ndims* nlandmarks];
l=l.(IdentityMatrix[ndims nlandmarks]-Mean[l].Mean[l]);
Return[l]]


(*  The CongruenceCoefficient[] function calculates a correlation matrix among a set of Procrustes superimposed 
  2D or 3D landmarks the sum of the covariances between the two multidimensional landmarks over all specimens 
  in the dataset divided by their pooled variance. Superimposed landmarks should be submitted in a three-
   dimensional array where the first dimension is the number of specimens, the second dimension is the number 
   of landmarks,and the third dimension is the number of coordinates in each landmark.  This function is used in 
   the Modularity[] function below. Note that data must be partitioned after using the Procrustes function above before being 
   entered in CongruenceCoefficient[] or CanonicalCorrelation[].  To partition data into 3D coordinates, use 
   Partition[Partition[Flatten[superimposed],ndims],nlandmarks] as detailed in the user guide.
 *)
CongruenceCoefficient[Superimposed_] := Module[ {i,j,N,stdevs},
N = Table[dotcvmentry[Superimposed,i,j], {i,1,Dimensions[Superimposed][[2]]}, {j,1,Dimensions[Superimposed][[2]]}];
stdevs=Sqrt[Tr[N,List]];
Return[N/Table[Table[stdevs[[i]]*stdevs[[j]],{j,Length[stdevs]}],{i,Length[stdevs]}]]
];

dotcvmentry[Superimposed_, col1_, col2_] := Module[ {i,n,s1,s2,p},
   n = 0;
If[Dimensions[Superimposed][[3]]==2,
   s1 =s2= {0,0},s1 =s2= {0,0,0}];

   For[i = 1, i <= Dimensions[Superimposed][[1]], i++,
      If[ !StringQ[Superimposed[[i,col1]]] && !StringQ[Superimposed[[i,col2]]],
         n++;
         s1 = s1 + Superimposed[[i,col1]];
         s2 = s2 + Superimposed[[i,col2]];
      ]
   ];
   If[n<=1,
      Print["there are too many missing data to covary columns ", col1, " and ", col2];
      Abort[];
   ];
   s1 = s1/n;
   s2 = s2/n;

   p = 0;
   For[i = 1, i <= Dimensions[Superimposed][[1]], i++,
      If[ !StringQ[Superimposed[[i,col1]]] && !StringQ[Superimposed[[i,col2]]],
         p = p + (Superimposed[[i,col1]] - s1) . (Superimposed[[i,col2]] - s2);
      ];
   ];
   p/(n-1)
];



(*  The CanonicalCorrelation[] function calculates a correlation matrix among a set of Procrustes superimposed 
  2D or 3D landmarks as the correlation between the major axis scores of the respective landmarks. Superimposed 
  landmarks should be submitted in a three-dimensional array where the first dimension is the number of 
  specimens, the second dimsension is the number of landmarks,and the third dimension is the number of coordinates 
   in each landmark. This function is used in the Modularity[] function below.
   
   	Updated 7 April 2018 to be consistent with current versions of Mathematica.  
 *)

CanonicalCorrelation[Superimposed_]:=Module[{pcs,i},
pcs=Transpose[Table[Transpose[PrincipalComponents[#-Mean[Superimposed[[1;;,i]]]&/@Superimposed[[1;;,i]]]][[1]],{i,Dimensions[Superimposed][[2]]}]];
Return[Correlation[pcs]];
];


(*  The Magwene[] function performs a graphical model analysis as described by Magwene 2001. "New tools for studying 
    integration and modularity", Evolution, 55: 1734-1745. Input consists of a symmetric correlation matrix, a single 
   integer for the sample size, and a vector of vertex labels. There are two optional variables, "Criterion" and 
   "Threshold",which specify which matrix to use for edge labels and the threshold value for drawing a graph edge.  By 
  default the criterion is the p-value for edge exclusion deviance and the threshold is to draw edges for those 
  whose probability is 0.95 or higher.  Other criteria are "EED" for the edge exclusion deviance value itself, "ES" 
   for edge strength, and "PC" for the partial correlation matrix.  The edge labels are given from the same matrix 
   specified by "Criterion".  The function returns a graphical model for integration and modularity,
    the partial R-square values for each variable with respect to all others, the matrix of partial correlations 
   among variables, the edge exclusion deviance matrix, the matrix of p-values for edge exclusion deviance, and the edge 
   strength matrix. 
*)

Magwene[P_,SampSize_,VertexLabels_,Criterion_:"EEDP",Threshold_:0.95]:=Module[{\[CapitalOmega],\[Rho],i,j,RSquare,EED, EEDP, ES,GraphEdges,EdgeData},
\[CapitalOmega] = Inverse[P];
RSquare = (Tr[\[CapitalOmega],List]-1)/Tr[\[CapitalOmega],List];
\[Rho] = Table[Table[-\[CapitalOmega][[i,j]]/(Sqrt[\[CapitalOmega][[i,i]]*\[CapitalOmega][[j,j]]]),{j,Length[P]}],{i,Length[P]}];
Do[\[Rho][[i,i]]=1,{i,Length[\[Rho]]}];
EED = Table[Table[(-SampSize)*Log[1-(\[Rho][[i,j]]^2)],{i,Length[\[Rho]]}],{j,Length[\[Rho]]}];
EEDP=Chop[Table[Table[CDF[ChiSquareDistribution[1],EED[[i,j]]],{i,Length[\[Rho]]}],{j,Length[\[Rho]]}]];
ES=Table[Table[(-0.5)*Log[1-(\[Rho][[i,j]]^2)],{i,Length[\[Rho]]}],{j,Length[\[Rho]]}];
EdgeData=Switch[Criterion,
"EEDP",EEDP,
"EED",EED,
"ES",ES,
"PC",\[Rho]
];
GraphEdges=Flatten[Table[Table[{VertexLabels[[i]]->VertexLabels[[j]],EdgeData[[i,j]]},{j,i+1,Length[VertexLabels]}],{i,Length[VertexLabels]-1}],1];
Return[{{GraphPlot[Select[GraphEdges,#[[2]]>Threshold&],VertexLabeling->True]},{"RSquare -> ",RSquare},{"Partial correlation matrix -> ",\[Rho]},{"Edge Exclusion Deviance Matrix ->",EED},{"Edge Exclusion Deviance P-values ->",EEDP},{"Edge Strength Matrix -> ",ES}}];
]



(*  The Modularity[] function does a simple analysis of modularity on geometric morphometric landmark data.  The function requires the landmarks 
(formatted in a matrix where each row contains the x, y (, z) coordinates of all landmarks for a single specimen), the number of landmarks, the 
number of landmark dimensions (2 or 3), and labels for each landmark.  The program returns results based on both the congruence coefficient and the 
canonical correlation (see Goswami and Polly, 2010).  For each of the two coefficients the eigenvalues of the correlation matrix are returned in a 
barchart along with their standard deviation (see Pavlicev, Cheverud and Wagner, 2009), as well as a Ward's linkage dendrogram showing modular clusters 
among the landmarks.  

	Updated 7 April 2018 to be consistent with current versions of Mathematica.  

*)
<<HierarchicalClustering`
Modularity[Landmarks_,NLands_,NDims_,Labels_,Randomization_:"False"]:=Module[{superimposed,residuals,P1,P2,EVStdev1,Plot1,EVStdev2,Plot2,
NinetyFivePercentCutoff1,NinetyFivePercentCutoff2,Tree1,Tree2,EV1,EV2,BootEV1,BootEV2,MaxMinBootEV1,MaxMinBootEV2,MeanBootEV1,MeanBootEV2,
Rows,Columns,dims,NumMods1,NumMods2,RelEVStdev1,RelEVStdev2},
superimposed=Procrustes[Landmarks,NLands,NDims];
superimposed=Partition[Partition[Flatten[superimposed],NDims],NLands];
{Rows,Columns,dims}=Dimensions[superimposed];
P1=CongruenceCoefficient[superimposed];
P2=CanonicalCorrelation[superimposed];
EV1=Chop[Eigenvalues[P1]];
EVStdev1=Sqrt[Plus@@((Eigenvalues[P1]-Mean[Eigenvalues[P1]])^2)/MatrixRank[P1]];
RelEVStdev1=Sqrt[Plus@@((Eigenvalues[P1]-Mean[Eigenvalues[P1]])^2)/MatrixRank[P1]]/Sqrt[MatrixRank[P1]];
EV2=Chop[Eigenvalues[P2]];
EVStdev2=Sqrt[Plus@@((Eigenvalues[P2]-Mean[Eigenvalues[P2]])^2)/MatrixRank[P2]];
RelEVStdev2=Sqrt[Plus@@((Eigenvalues[P2]-Mean[Eigenvalues[P2]])^2)/MatrixRank[P2]]/Sqrt[MatrixRank[P2]];


If[Randomization=="RandomizationTest",

BootEV1=Table[Eigenvalues[CongruenceCoefficient[Partition[superimposed[[#[[1]],#[[2]]]]&/@Flatten[ Table[Table[{RandomInteger[{1,Rows}],RandomInteger[{1,Columns}]},{Columns}],{Rows}],1],NLands]]],{100}];
MaxMinBootEV1=Chop[Transpose[{Max[#]&/@Transpose[BootEV1],Min[#]&/@Transpose[BootEV1]}]];
MeanBootEV1=Chop[Mean[#]&/@Transpose[BootEV1]];
NumMods1=Length[Select[EV1-(#[[1]]&/@MaxMinBootEV1),Positive]];
Plot1=Graphics[{
{Gray,Table[Rectangle[{x-.45,0},{x+.45,EV1[[x]]}],{x,Length[EV1]}]},
{Black,AbsoluteThickness[2],Dashed,Table[Line[{{x,MaxMinBootEV1[[x,2]]},{x,MaxMinBootEV1[[x,1]]}}],{x,Length[MaxMinBootEV1]}]},
{Black,AbsolutePointSize[5],Table[Point[{x,MeanBootEV1[[x]]}],{x,Length[MeanBootEV1]}]}
},PlotLabel->"Eigenvalues"];
BootEV2=Table[Eigenvalues[CanonicalCorrelation[Partition[superimposed[[#[[1]],#[[2]]]]&/@Flatten[ Table[Table[{RandomInteger[{1,Rows}],RandomInteger[{1,Columns}]},{Columns}],{Rows}],1],NLands]]],{100}];
MaxMinBootEV2=Chop[Transpose[{Max[#]&/@Transpose[BootEV2],Min[#]&/@Transpose[BootEV2]}]];
MeanBootEV2=Chop[Mean[#]&/@Transpose[BootEV2]];
NumMods2=Length[Select[EV2-(#[[1]]&/@MaxMinBootEV2),Positive]];
Plot2=Graphics[{
{Gray,Table[Rectangle[{x-.45,0},{x+.45,EV2[[x]]}],{x,Length[EV2]}]},
{Black,AbsoluteThickness[2],Dashed,Table[Line[{{x,MaxMinBootEV2[[x,2]]},{x,MaxMinBootEV2[[x,1]]}}],{x,Length[MaxMinBootEV2]}]},
{Black,AbsolutePointSize[5],Table[Point[{x,MeanBootEV2[[x]]}],{x,Length[MeanBootEV2]}]}
},PlotLabel->"Eigenvalues"];
Tree1=DendrogramPlot[DirectAgglomerate[Chop[1-Abs[P1]],Labels,Linkage->"Ward"],LeafLabels->(#&),HighlightLevel->If[NumMods1==0,None,NumMods1]];
Tree2=DendrogramPlot[DirectAgglomerate[Chop[1-Abs[P2]],Labels,Linkage->"Ward"],LeafLabels->(#&),HighlightLevel->If[NumMods2==0,None,NumMods2]];
Return[Style[Grid[{
{Style["Modularity Results",FontWeight->Bold],""},
{Style["Congruence Coefficient",FontSlant->Italic],""},
{"Eignvalue SD: "<>ToString[EVStdev1]},
{"Relative Eigenvalue SD: "<>ToString[RelEVStdev1]},
{"Expectation of Random Rel. Eigenvalue SD: "<>ToString[Sqrt[1/Length[Landmarks]//N]]},
{"Estimated Number of Modules: "<>ToString[NumMods1]},
{Plot1,Tree1},
{Style["Canonical Correlation",FontSlant->Italic],""},
{"Eignvalue SD: "<>ToString[EVStdev2]},
{"Relative Eigenvalue SD: "<>ToString[RelEVStdev2]},
{"Expectation of Random Rel. Eigenvalue SD: "<>ToString[Sqrt[1/Length[Landmarks]//N]]},
{"Estimated Number of Modules: "<>ToString[NumMods2]},
{Plot2,Tree2}
},Alignment->Left],FontFamily->"Helvetica"]
]
,

Plot1=BarChart[EV1,PlotLabel->"Eigenvalues"]; Plot2=BarChart[EV2,PlotLabel->"Eigvenvalues"];
Tree1=DendrogramPlot[DirectAgglomerate[Chop[1-Abs[P1]],Labels,Linkage->"Ward"],LeafLabels->(#&)];
Tree2=DendrogramPlot[DirectAgglomerate[Chop[1-Abs[P2]],Labels,Linkage->"Ward"],LeafLabels->(#&)];
Return[Style[Grid[{
{Style["Modularity Results",FontWeight->Bold],""},
{Style["Congruence Coefficient",FontSlant->Italic],""},
{"Eigenvalue SD: "<>ToString[EVStdev1]},
{"Relative Eigenvalue SD: "<>ToString[RelEVStdev1]},
{"Expectation of Random Rel. Eigenvalue SD: "<>ToString[Sqrt[1/Length[Landmarks]//N]]},
{Plot1,Tree1},
{Style["Canonical Correlation",FontSlant->Italic],""},
{"Eigenvalue SD: "<>ToString[EVStdev2]},
{"Relative Eigenvalue SD: "<>ToString[RelEVStdev2]},
{"Expectation of Random Rel. Eigenvalue SD: "<>ToString[Sqrt[1/Length[Landmarks]//N]]},
{Plot2,Tree2}
},Alignment->Left],FontFamily->"Helvetica"]
]

];



]



(*
 * Usage:  mantel[A,B,n], where A and B are square matrices of the
 * same size, and n is an integer.  Returns the matrix correlation 
 * between the two original matrices, the fraction of times in the
 * n Monte Carlo randomization runs that the matrix correlation is 
 * higher than a random correlation, and the probability that the 
 * correlation is greater than random.
 *)

randperm[n_] := Module[ {l,t,r},
    l = Table[i, {i,1,n}];
    For [i = n, i > 1, i--,
	r = Random[Integer, {1,i}];
	t = l[[i]];
	l[[i]] = l[[r]];
	l[[r]] = t;
    ];
    Return[l];
];

randmperm[A_] := Module[ {n,p},
    If [!MatrixQ[A, NumberQ],
	Print["fatal error in mantel: bad matrix"];
	Abort[];
    ];
    If [Dimensions[A][[1]] != Dimensions[A][[2]],
	Print["fatal error in mantel: matrix not square"];
	Abort[];
    ];

    n = Dimensions[A][[1]];
    p = randperm[n];
    Return[Table[ A[[p[[i]],p[[j]]]], {i,1,n}, {j,1,n} ]];
];

mantelcor[l1_,l2_] := (l1.l2 - (Plus @@ l1) * (Plus @@ l2) / Length[l1]) /
(Sqrt[l1.l1 - (Plus @@ l1) * (Plus @@ l1) / Length[l1]] *
 Sqrt[l2.l2 - (Plus @@ l2) * (Plus @@ l2) / Length[l2]]);

mantelmcor[A_,B_] := mantelcor[Flatten[A],Flatten[B]];

Mantel[A_,B_,n_] := Module[ {i,t,tt,npass},
	t = mantelmcor[A,B];
	npass = 0.0; 
	For [i = 1, i <= n, i++,
		tt = mantelmcor[A,randmperm[B]];
		If [t >=  tt, npass++]
	 ];
    Return[{t, npass/n, 1-npass/n}];
];



(*Subsample[A,m]- subsamples data matrix A down to n rows*)
RP[m_]:=Module[{rpx,rpi,rpj,rpt},
rpx=Table[rpi,{rpi,1,m}];
For[rpi=1,rpi<=m,rpi=rpi+1,rpj=Random[Integer,{1,m}];
rpt=rpx[[rpi]];
rpx[[rpi]]=rpx[[rpj]];
rpx[[rpj]]=rpt;];
Return[rpx];
];

Subsample[A_,m_]:=Module[{rarx,rari,rary},
rarx=RP[Length[A]];
rary=Table[rarx[[rari]],{rari,1,m}];
Return[A[[rary]]]
];



(*subsamples matrix A from n-1 to n-m*)
(*computes correlation coefficient with constant matrix B, can be same or different matrix*)
(*k=number of landmarks, d=number of dimensions*)
(*m=max.number of rows to remove*)
(*t=no.of trials at each subsampling step*)

Diag[P_]:=Module[{i,j},Flatten[Table[Table[P[[i,j]],{j,i+1,Length[P]}],{i,Length[P]-1}]]]

SubsampleMatrix[A_,B_,nland_,ndims_,m_,t_]:=Module[{i,j,l},
l=Length[A];
Return[Table[
Correlation[Diag[CongruenceCoefficient[Partition[Partition[Flatten[
Procrustes[Subsample[A,l-i],nland,ndims]],ndims],nland]]],Diag[CongruenceCoefficient[Partition[Partition[Flatten[
Procrustes[B,nland,ndims]],ndims],nland]]]],
{i,1,m},{j,1,t}]]
];
