/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
   clusterMHDFoam

Description
    Density-based compressible MHD flow solver


Author: Carlos Xisto


All rights reserved.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "string.H"
#include "OFstream.H"


//Numerical flux functions for MHD
#include "AUSMpw-MHD.H" //Han2009
#include "AUSMup-MHD.H" //Liou 2006 version modified for MHD
#include "TVDLF-MHD.H"  //Thoth96
#include "Kurganov-MHD.H"
#include "HLLC-L.H"     //Li2005
#include "HLLD.H"       //Miyoshi2005
//#include "Roe-MHD.H"


//Numerical flux functions for Gas dynamics and Low Mag-Re formulation
#include "HLLC.H"
#include "Roe.H"
#include "AUSMupPlus.H"
#include "Kurganov.H"
#include "TVDLF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
  #include "setRootCase.H"

  #include "createTime.H"
  #include "createMesh.H"

  #include "fluxMHD.H"

  #include "createFields.H"
  #include "readTimeControls.H"

  #include "CourantNo.H"
  #include "setInitialDeltaT.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




  word RK4("yes");
  if (mesh.schemesDict().readIfPresent("RK4", RK4))
  {
    if
    (
      (RK4 == "yes") ||
      (RK4 == "no")

    )

    {
      Info<< "4-Stage RK on?  " << RK4 << nl
      << endl;
    }

  }



  Info<< "\nStarting time loop\n" << endl;
  while (runTime.run())
  {



    // Compute the numerical flux across all faces
    #include "numericalFlux.H"


    #include "CourantNo.H"
    #include "readTimeControls.H"
    #include "setDeltaT.H"

    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;



    volScalarField rho0=rho;
    volVectorField rhoU0=rhoU;
    volVectorField B0=B;
    volScalarField pB0=pB;
    volScalarField rhoE0=rhoE;

    if (RK4 == "yes")
    {
      /*
      =============================================================
      This is a low-storage second-order version of the R-K method
      Only the last residual is stored

      U(0) = U(n)
      U(1) = U(0) - alpha(1)*R[U(0)]
      ...
      ...
      U(k) = U(0) - alpha(k)*R[U(k-1)]
      U(n+1) = U(k)

      alpha(k) = 1
      alpha(k-1) = 1/2 (Required for second order)

      Note: If only one stage is selected (alpha(1,2,3)=0)
      the Euler explicit scheme is recovered
      -------------------------------------------------------------
      see the Book of Blazek pp. 182. for more coefs
      -------------------------------------------------------------
      =============================================================
      */


      //4Stage Runge-Kutta coefs
      scalarList RK4alpha(0);
      RK4alpha.append(0.125);
      RK4alpha.append(0.25);
      RK4alpha.append(0.5);


      //Start Runge-Kutta loop for calculating the fluxes
      for (int i = 0; i < RK4alpha.size(); i++ )
      {
        //Compute the numerical flux across all faces
        #include "numericalFlux.H"

        Info<< "RKStage:: " << i+1 << endl;

        //solve(fvm::ddt(rho) == -RK4alpha[i]*fvc::div(phi));

        rho  = rho0  - RK4alpha[i]*runTime.deltaT()*fvc::div(phi);

        if(magnetic)
        {
          #include "BEqnRK4.H"    //Induction equation
        }
        #include "UEqnRK4.H"    //momentum equation
        #include "eEqnRK4.H"    //Energy equation
      }

    }


    #include "numericalFlux.H"
    //initializing the variables for computing the last RK step (alpha(k)=1)
    rho=rho0;
    rhoU=rhoU0;
    B=B0;
    pB=pB0;
    rhoE=rhoE0;

    solve(fvm::ddt(rho) + fvc::div(phi)); //Solve continuity

    if(magnetic)
    {
      #include "BEqn.H"    //Induction equation
    }

    #include "UEqn.H"    //momentum equation
    #include "eEqn.H"    //Energy equation


    runTime.write();

    divBf = fvc::div(Bf); //Update div(Bf)

    //Number of cells for average (works in parallel)
    label Ncells = returnReduce(mesh.cells().size(), sumOp<label>());

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl
    << "max(div(Bf)) = " << gMax(mag(divBf.internalField()))
    << "  average(div(Bf)) = " << gSum(mag(divBf.internalField()))/Ncells
    << nl
    << nl << endl;
  }

  Info<< "End\n" << endl;
  return 0;
}
