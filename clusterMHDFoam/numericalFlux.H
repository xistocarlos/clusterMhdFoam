
{


  /*____________________________________________________________________________
  LEFT & RIGTH FACE VALUES

  Using J. Greenshields approach

  Christopher J. Greenshields, H- G. Weller, Luca Gasparini and J. M. Reese
  "Implementation of semi-discrete, non-staggered central schemes in a colocated,
  polyhedral, finite volume framework,for high-speed viscous flows", IJNMF,2009.
  ______________________________________________________________________________
  */
  surfaceVectorField    rhoUL = fvc::interpolate(rhoU, pos, "reconstruct(U)");
  surfaceVectorField    rhoUR = fvc::interpolate(rhoU, neg, "reconstruct(U)");

  surfaceVectorField    BL = fvc::interpolate(B, pos, "reconstruct(B)");
  surfaceVectorField    BR = fvc::interpolate(B, neg, "reconstruct(B)");

  surfaceScalarField    rhoL = fvc::interpolate(rho, pos, "reconstruct(rho)");
  surfaceScalarField    rhoR = fvc::interpolate(rho, neg, "reconstruct(rho)");

  surfaceScalarField    pL = fvc::interpolate(p, pos, "reconstruct(p)");
  surfaceScalarField    pR = fvc::interpolate(p, neg, "reconstruct(p)");

  surfaceScalarField    rhoEL = fvc::interpolate(rhoE, pos, "reconstruct(T)");
  surfaceScalarField    rhoER = fvc::interpolate(rhoE, neg, "reconstruct(T)");

  surfaceScalarField    pBL = fvc::interpolate(pB, pos, "reconstruct(pB)");
  surfaceScalarField    pBR = fvc::interpolate(pB, neg, "reconstruct(pB)");

  surfaceScalarField gammaL = fvc::interpolate(gamma, pos, "reconstruct(gamma)");
  surfaceScalarField gammaR = fvc::interpolate(gamma, neg, "reconstruct(gamma)");

  surfaceVectorField    UL = rhoUL/rhoL;
  surfaceVectorField    UR = rhoUR/rhoR;


  if(fluxScheme == "AUSMpw-MHD")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      AUSMpwFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        pBL[f],   pBR[f],
        BL[f],    BR[f],
        Bf[f], phi[f],
        phiUpB[f], phiUB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        AUSMpwFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          pBL.boundaryField()[iPatch][f],
          pBR.boundaryField()[iPatch][f],
          BL.boundaryField()[iPatch][f],
          BR.boundaryField()[iPatch][f],
          Bf.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiUB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }


  else if (fluxScheme == "Kurganov-MHD")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      KurganovMHDFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        pBL[f],   pBR[f],
        BL[f],    BR[f],
        Bf[f], phi[f],
        phiUpB[f], phiUB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        KurganovMHDFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          pBL.boundaryField()[iPatch][f],
          pBR.boundaryField()[iPatch][f],
          BL.boundaryField()[iPatch][f],
          BR.boundaryField()[iPatch][f],
          Bf.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiUB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }

  else  if(fluxScheme == "AUSMup-MHD")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      AUSMupMHDFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f], rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],   UR[f],
        pL[f],   pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f], Ms[f],
        pBL[f],   pBR[f],
        BL[f],   BR[f],
        Bf[f], phi[f],  phiUpB[f],
        phiUB[f], phiEpB[f], cfmax[f]
      );
    }


    // Boundary faces
    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")


      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {
        AUSMupMHDFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          Ms.boundaryField()[iPatch][f],
          pBL.boundaryField()[iPatch][f],
          pBR.boundaryField()[iPatch][f],
          BL.boundaryField()[iPatch][f],
          BR.boundaryField()[iPatch][f],
          Bf.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiUB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }





  else if (fluxScheme == "TVDLF-MHD")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      TVDLFMHDFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        pBL[f],   pBR[f],
        BL[f],    BR[f],
        Bf[f], phi[f],
        phiUpB[f], phiUB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        TVDLFMHDFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          pBL.boundaryField()[iPatch][f],
          pBR.boundaryField()[iPatch][f],
          BL.boundaryField()[iPatch][f],
          BR.boundaryField()[iPatch][f],
          Bf.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiUB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }

  else if (fluxScheme == "HLLC-L")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      HLLCLFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        pBL[f],   pBR[f],
        BL[f],    BR[f],
        Bf[f], phi[f],
        phiUpB[f], phiUB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        HLLCLFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          pBL.boundaryField()[iPatch][f],
          pBR.boundaryField()[iPatch][f],
          BL.boundaryField()[iPatch][f],
          BR.boundaryField()[iPatch][f],
          Bf.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiUB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }


  else if (fluxScheme == "HLLD")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      HLLDFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        pBL[f],   pBR[f],
        BL[f],    BR[f],
        Bf[f], phi[f],
        phiUpB[f], phiUB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        HLLDFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          pBL.boundaryField()[iPatch][f],
          pBR.boundaryField()[iPatch][f],
          BL.boundaryField()[iPatch][f],
          BR.boundaryField()[iPatch][f],
          Bf.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiUB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }

  /*
  ==================================================================================
  Flux functions for MHD - End
  ==================================================================================
  */

  //The following fuctions are only working with gas-dynamic or low-Mag Re
  //formulations.
  





  if (fluxScheme == "HLLC")
  {
    Info<< "Warning: This is a Gas-Dynamic only scheme " << endl;
    // Interior faces
    forAll (mesh.Sf(), f)
    {
      HLLCFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        phi[f],  phiUpB[f],
        phiEpB[f], cfmax[f]
      );
    }


    // Boundary faces
    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")


      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {
        HLLCFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }

    }

  }
  else if (fluxScheme == "Roe")
  {

    Info<< "Warning: This is a Gas-Dynamic only scheme " << endl;
    // Interior faces
    forAll (mesh.Sf(), f)
    {
      RoeFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        phi[f],  phiUpB[f],
        phiEpB[f], cfmax[f]
      );
    }


    // Boundary faces
    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")


      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {
        RoeFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }


  else  if(fluxScheme == "AUSMupPlus")
  {

    Info<< "Warning: This is a Gas-Dynamic only scheme " << endl;

    // Interior faces
    forAll (mesh.Sf(), f)
    {
      AUSMupPlusFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f], rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],   UR[f],
        pL[f],   pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f], Ms[f],
        phi[f],  phiUpB[f],
        phiEpB[f], cfmax[f]
      );
    }


    // Boundary faces
    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")


      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {
        AUSMupPlusFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          Ms.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }


  else if (fluxScheme == "Kurganov")
  {
    Info<< "Warning: This is a Gas-Dynamic only scheme " << endl;

    // Interior faces
    forAll (mesh.Sf(), f)
    {
      KurganovFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        phi[f],   phiUpB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        KurganovFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }


  else if (fluxScheme == "TVDLF")
  {


    // Interior faces
    forAll (mesh.Sf(), f)
    {
      TVDLFFlux (
        mesh.Sf()[f],
        mesh.magSf()[f],
        rhoL[f],  rhoR[f],
        rhoUL[f], rhoUR[f],
        UL[f],    UR[f],
        pL[f],    pR[f],
        rhoEL[f], rhoER[f],
        gammaL[f], gammaR[f],
        phi[f],  phiUpB[f],
        phiEpB[f], cfmax[f]
      );
    }



    forAll (mesh.boundaryMesh(), iPatch)
    {
      if(mesh.boundaryMesh().physicalTypes()[iPatch] != "empty")
      forAll (mesh.Sf().boundaryField()[iPatch], f)
      {

        TVDLFFlux (
          mesh.Sf().boundaryField()[iPatch][f],
          mesh.magSf().boundaryField()[iPatch][f],
          rhoL.boundaryField()[iPatch][f],
          rhoR.boundaryField()[iPatch][f],
          rhoUL.boundaryField()[iPatch][f],
          rhoUR.boundaryField()[iPatch][f],
          UL.boundaryField()[iPatch][f],
          UR.boundaryField()[iPatch][f],
          pL.boundaryField()[iPatch][f],
          pR.boundaryField()[iPatch][f],
          rhoEL.boundaryField()[iPatch][f],
          rhoER.boundaryField()[iPatch][f],
          gammaL.boundaryField()[iPatch][f],
          gammaR.boundaryField()[iPatch][f],
          phi.boundaryField()[iPatch][f],
          phiUpB.boundaryField()[iPatch][f],
          phiEpB.boundaryField()[iPatch][f],
          cfmax.boundaryField()[iPatch][f]
        );
      }
    }
  }
}
