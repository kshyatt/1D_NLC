/*******************************************************************
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
baesd on a Lanczos code from November 2007
********************************************************************/

#include <utility>
#include <fstream> 
#include <vector> 
#include <math.h>
using namespace std;

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <string>

#include "CPU/Lanczos_07.h"
#include "CPU/GenHam.h"
#include "CPU/simparam.h"
#include "CPU/magnetization.h"
#include "graphs.h"
#include <mpi.h>

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int CurrentArg = 1;
    string InputFile;
    string OutputFile = "output_2d.dat";
    while (CurrentArg < argc)
    {
        if (argv[ CurrentArg ] == string("-i") || argv[ CurrentArg ] == string("--input") )
        {
            InputFile = string(argv[ CurrentArg + 1 ]);
        }
        if (argv[ CurrentArg ] == string("-o") || argv[ CurrentArg ] == string("--output"))
        {
            OutputFile = string(argv[ CurrentArg + 1 ]);
        }
        CurrentArg++;
    }
    
    double energy;
    double chi;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J = 1.;

    vector <long double> eVec;

    vector< graph > fileGraphs; //graph objects
    
    vector<double> EnergyWeightHigh;
    vector<double> MagnetizationWeightHigh;
    
    ReadGraphsFromFile(fileGraphs, InputFile);
    
    if ( rank == 0 ) 
    {
        ofstream fout(OutputFile.c_str());
        fout.precision(10);
        double* EnergyResults = (double*)malloc((size - 1)*sizeof(double));
        double* MagnetizationResults = (double*)malloc((size - 1)*sizeof(double));
        for( int i = 0; i < size - 1; i++ )
        {
            MPI_Recv(EnergyResults + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &status); //grab results from other processes
            MPI_Recv(MagnetizationResults + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &status); //grab results from other processes
            fout<<"h = "<<(i + 1)/(10*J) + 1.<<" J = "<<J<<" Energy = "<<EnergyResults[i]<<" Magnetization = "<<MagnetizationResults[i]<<endl;
        }
        fout.close();
    }
    
    double h = (double)rank/(10*J) + 1.;
      
    EnergyWeightHigh.push_back(-h); //Weight for site zero
    double EnergyRunningSumHigh = WeightHigh[0];      
    MagnetizationWeightHigh.push_back(1.); //Weight for site zero
    double MagnetizationRunningSumHigh = MagnetizationWeightHigh[0];      
    
    if ( rank )
    {
        for (int i = 1; i < fileGraphs.size(); i++)
        { //skip the zeroth graph
	
	    //---High-Field---
	        GENHAM HV(fileGraphs.at(i).NumberSites, J, h, fileGraphs.at(i).AdjacencyList, fileGraphs.at(i).LowField); 

            LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
    
            HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos
            energy = lancz.Diag(HV, 1, 1, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
            chi = Magnetization(eVec, fileGraphs.at(i).NumberSites);
            EnergyWeightHigh.push_back(energy);
            MagnetizationWeightHigh.push_back(chi);
    
            for (int j = 0; j < fileGraphs[ i ].SubgraphList.size(); j++)
            {
	            EnergyWeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * EnergyWeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];
	            MagnetizationWeightHigh.back() -= fileGraphs[ i ].SubgraphList[ j ].second * MagnetizationWeightHigh[ fileGraphs[ i ].SubgraphList[ j ].first ];
            }

            EnergyRunningSumHigh += fileGraphs[ i ].LatticeConstant * EnergyWeightHigh.back();
            MagnetizationRunningSumHigh += fileGraphs[ i ].LatticeConstant * MagnetizationWeightHigh.back();
	
        }
        MPI_Send( &EnergyRunningSumHigh, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //send resulting sum to controller process
        MPI_Send( &MagnetizationRunningSumHigh, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //send resulting sum to controller process

    }

    return 0;

}
