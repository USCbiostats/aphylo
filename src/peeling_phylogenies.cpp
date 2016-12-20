#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <random4f.h>
#include <cholesky.cxx>
#include <cholesky21.h>

const	int		N		=	326,	// number of nodes
				L		=	175,	// number of leaves
				F		=	3,		// number of gene functions
				FF		=	8,		// number of possible joint function states (2^F)
				MaxOff	=	10,		// maximum number of offspring nodes for any one internal node
				P	=	5;
				
const	double	qinit[2]	=	{0.020,0.010},		// initial misclassification probabilities (FP,FN)
				P0init		=	0.1,				// initial functional status of root node
				muinit[2]	=	{0.004,.001};		// initial mutation rates (loss,gain)



int	f,l,n,o,p,p1,p2,par,xx,
	parent[N],			// parent node of node n
	Noff[N],			// number of offspring node from node n
	OffId[N][MaxOff],	// IDs of offspring of node n
	NodeId[N],			// input array of node IDs (original numbering system) 
	ParentId[N],		// input array of parent node IDs (ditto)
	NodeType[N],		// 0 = duplication, 1 = speciation, 2 = leaf
	Z[N][F],			// experimentally determined function data for leaf nodes (9 = missing)
	x[FF];				// joint function state vector

double	prX[N][FF],			// peeled node probabilities based on all decendent nodes (including leaves)
		psi[2],mu[2],pi[2],		// current values of parameters
		parm[P],				// vector of all psi,mu,pi parameters 
		U[P],					// score vector
		Info[P][P],InvInfo[P][P];	// information matrix and its inverse
		
float	BranchLength[N];		// input array of branch lengths for node n
		

FILE *sum,*dag,*dat,*itr;


void Initialize()
{	psi[0] = qinit[0]; psi[1] = qinit[1];
	mu[0] = muinit[0]; mu[1] = muinit[1]; 
	pi[0] = 1-P0init;  pi[1] = P0init;

}

void InputData()
{	dag = fopen("PTHR11848.dag","r");
	dat = fopen("PTHR11848 sorted.txt","r");
	itr = fopen("phylogene.itr","w");

	for (n=0; n<N; n++)	fscanf (dag,"%d %d %d",&NodeId[n],&NodeType[n],&ParentId[n]);;

	// store experimental data for leaf nodes
	//		note: this file includes some IDs that do not appear in the .dag file
	//				so these is simply skipped

	n=0; int LeafId,ZZ[F];
	while (n<L)
	{	for (f=0; f<F; f++) fscanf (dat,"%d",&ZZ[f]);	
		fscanf (dat,"%d",&LeafId);
		while (NodeId[n]<LeafId) 
		{	for (f=0; f<F; f++) Z[n][f] = 8;	// no experimental data for internal nodes
			n++;
		}
		if (NodeId[n]==LeafId)
		{	for (f=0; f<F; f++) Z[n][f] = ZZ[f];
			n ++;
		}
	}

	// reassign IDs and parent pointers in consequetive order
	parent[0] = -1;	// root node
	for (n=1; n<N; n++)
	{	int found=0,par=0;
		while (!found)
		{	if (ParentId[n] == NodeId[par]) 
			{	parent[n] = par;
				found = 1;
			}
			par ++;
		}
	}
}


void GetPeelingSequence()
{	int nn=0;
	memset(Noff,0,sizeof(Noff));
	n=0;
	while (n<N)
	{	
		// first check to see if the parent node has already been identified.
		//	if so, add it to the list of that parent's offspring 
		
		int found=0;
		for (par=0; par<n; par++)
			if (parent[n]==par) 
			{	found=1;
				OffId[par][Noff[par]] = n;
				Noff[par] ++;
			}

		// otherwise it has a new parent, so set up a pointer to it as the first offspring of that one

		if (!found && n)
		{	par = parent[n];
			Noff[par] = 1;
			OffId[par][0] = n;
		}
//			for (o=n+1; o<N; o++)
//				if (parent[o]==p)
//				{	Noff[par] ++;
//					OffId[par][Noff[par]-1] = o;
//				}
//		}
		n ++;
	}
}

void GetX(int xx)
{	// decomposes index of joint function states xx into vector of single-function states x[FF]
	for (f=0; f<F; f++)
	{	x[f] = xx%2;
		xx /= 2;
	}	
}


void MutationProbabilities(int nodetype, double time)
{
	// this one still needs to be developed
	//		to compute mutation probabilities as a function of node type and branch length
	//		and also for vectors of correlated functions
}

double LogLike(double parm[P])
{	
	// copy parameters from the generic parameter vector parm into psi, mu, pi arrays
	
	psi[0] = parm[0];	psi[1] = parm[1];
	mu[0]  = parm[2];	mu[1]  = parm[3]; 
	pi[0]  = 1-parm[4];	pi[1]  = parm[4];

	// leaf node probabilities: 
	//		if function has been experimentally determined, 
	//			then use the misclassification probabilities
	//		otherwise 1 

	for (n=0; n<N; n++)
		if (Noff[n]==0)			// perform this step only for leaf nodes
			for (int xxleaf=0; xxleaf<FF; xxleaf++)		// loop over joint function states
			{	GetX(xxleaf);		// decompose xx into vector of function-specific states x[FF]
				prX[n][xxleaf] = 1;
				for (f=0; f<F; f++)
				{	if (Z[n][f] != 9) 	// include only functions with experimental data
					{	if (Z[n][f] == x[f]) 	prX[n][xxleaf] *= 1-psi[x[f]];
										else	prX[n][xxleaf] *=   psi[x[f]];
					}
				}
			}
	
	
	// internal node probabilities
	//		if gain of function, use mu[0]
	//		else if loss of function, use mu[1]
	//		else if no change, use 1-mu[0] or 1-mu[1] depending on parental state
	//	note: for now, I'm treating mutations as independent across functions
	//			with probabilities that will eventually depend on node type and branch lengths.
	//		For now, these probabilities are stored in the mu array,
	//			but eventually they will be functions of other parameters 
	//			that will be stored in the parm array
	//	In this file, nodes are arranged with the root first, then its descendents
	//		so the parents of each node appear before their offspring in the file
	//		Hence, the peeling must start at the bottom of the file and work upwards towards the root
	
	for (n=N-1; n>=0; n--) 
		if (Noff[n])
		{	for (int xxpar=0; xxpar<FF; xxpar++)
			{	prX[n][xxpar] = 1;
				int xpar[F]; GetX(xxpar); for (f=0; f<F; f++) xpar[f] = x[f];
				for (o=0; o<Noff[n]; o++)
				{	int offid = OffId[n][o];
					MutationProbabilities(NodeType[n],BranchLength[offid]);
					double sumprob=0;
					for (int xxoff=0; xxoff<FF; xxoff++)
					{	int xoff[F]; GetX(xxoff); 
						double prmut=1;
						for (f=0; f<F; f++) 
						{	xoff[f] = x[f];
							if (xpar[f]==0)					// parent node doesn't have function
							{	if (xoff[f]==1)	prmut *=   mu[0];	// gain of function mutation
										else	prmut *= 1-mu[0];	// no gain
							}
									else						// parent node does have function
							{	if (xoff[f]==0) prmut *=   mu[1];	// loss of function mutation
										else	prmut *= 1-mu[1];  // no loss
							}
						}
						sumprob += prmut * prX[offid][xxoff];
					}
					prX[n][xxpar] *= sumprob;
				}
			}	
		}
	
	// root node probabilities

	double like=1;
	for (int xxroot=0; xxroot<FF; xxroot++)
	{	GetX(xxroot); 
		double prFF0=1; for (f=0; f<F; f++) prFF0 *= pi[x[f]];
		like += prFF0 * prX[0][xxroot];
	}
	double loglike = log(like);
	return (loglike);
}


double NumericalDerivatives()
{	double dparm = 0.0001;
	double	parmP0[P],parmM0[P],
			parmPP[P],parmPM[P],parmMP[P],parmMM[P];
	for (p1=0; p1<P; p1++)
	{	parmP0[p1] = parm[p1]; parmM0[p1] = parm[p1];
		parmPP[p1] = parm[p1]; parmPM[p1] = parm[p1];
		parmMP[p1] = parm[p1]; parmMM[p1] = parm[p1];
	}
	double LL0 = LogLike(parm); 
	for (p1=0; p1<P; p1++)
	{	parmP0[p1] += dparm; parmM0[p1] -= dparm;
		double LLP = LogLike(parmP0); double LLM = LogLike(parmM0); 
		U[p1] = (LLP - LLM)/(2*dparm);
		Info[p1][p1] = (2*LL0 - LLP - LLM) / (dparm*dparm);
		parmP0[p1] -= dparm; parmM0[p1] += dparm;
		for (p2=0; p2<P; p2++)
			if (p1 != p2)
			{	parmPP[p1] += dparm; parmPP[p2] += dparm;
				parmPM[p1] += dparm; parmPM[p2] -= dparm;
				parmMP[p1] -= dparm; parmMP[p2] += dparm;
				parmMM[p1] -= dparm; parmMM[p2] -= dparm;
				double 	LLPP = LogLike(parmPP); double LLPM = LogLike(parmPM); 
				double	LLMP = LogLike(parmMP); double LLMM = LogLike(parmMM); 
				Info[p1][p2] = (LLPP - LLPM - LLMP + LLMM) / (4*dparm*dparm);
				parmPP[p1] -= dparm; parmPP[p2] -= dparm;
				parmPM[p1] -= dparm; parmPM[p2] += dparm;
				parmMP[p1] += dparm; parmMP[p2] -= dparm;
				parmMM[p1] += dparm; parmMM[p2] += dparm;
			}
	}
	return (LL0);
}


void MaximizeLikelihood()
{	// store specific parameters psi, mu, pi into generic array parm
	
	parm[0] = psi[0]; parm[1] = psi[1];
	parm[2] = mu[0]; parm[3] = mu[1];
	parm[4] = pi[1];
	
	int conv=0; double chisq;
	while (!conv)
	{	double LL = NumericalDerivatives();
		int error = InvertPDS(Info[0],P,InvInfo[0]);
		if (!error)
		{	chisq=0;
			for (p1=0; p1<P; p1++)
			for (p2=0; p2<P; p2++)
			{	parm[p1] += InvInfo[p1][p2] * U[p2];
				chisq += U[p1] * InvInfo[p1][p2] * U[p2];
			}
			if (chisq < 0.001) conv=1;
		}
		else
		{	printf ("\nERROR: non-positive-definite information matrix ... terminating");
			for (p1=0; p1<P; p1++)
			{	printf ("\n%d  %8.4f  ",p1,Info[p1][p1]);
				for (p2=0; p2<P; p2++)
					printf (" %6.3f",Info[p1][p2]/sqrt(Info[p1][p1]*Info[p2][p2]));
			}
			conv=1;
			chisq=0;
		}
		printf ("\n%6.3f %8.3f ",chisq,LL);
		for (p=0; p<P; p++) printf (" %7.4f",parm[p]);
		fprintf (sum,"\n%6.3f %8.3f ",chisq,LL);
		for (p=0; p<P; p++) fprintf (sum," %7.4f",parm[p]);
	}
	psi[0] = parm[0];	psi[1] = parm[1];
	mu[0]  = parm[2];	mu[1]  = parm[3]; 
	pi[0]  = 1-parm[4];	pi[1]  = parm[4];
}

void Analyze()
{	GetPeelingSequence();
	int conv=0;
	while (!conv)
	{	MaximizeLikelihood();
	}

}

void Tabulate()
{

}

void Summarize()
{	
	fprintf (sum,"\n\nSIMULATION SUMMARY\nstat  mean   (SD)");
}

void main()
{	sum = fopen("Peeling phylogenies.sum","w");
	Initialize();
	InputData();
	Analyze();
	Summarize();
	fclose(sum);
}
