//==========================================================
dimensionedScalar Minf = 0.0;

surfaceScalarField Ms
(
  IOobject
  (
    "Ms",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("Ms", Minf)
);



if (fluxScheme == "AUSMup-MHD")
{
  IOdictionary AUSM
  (
    IOobject
    (
      "AUSM",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  Info<< "AUSM file read, check Minf!\n" << nl << endl;

  Minf=AUSM.lookup("Minf");
  Info<< "Minf = " << Minf.value() <<  nl << endl;
  Ms = Minf;
}

else if (fluxScheme == "AUSMupPlus")
{
  IOdictionary AUSM
  (
    IOobject
    (
      "AUSM",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  );

  Info<< "AUSM file read, check Minf!\n" << nl << endl;

  Minf=AUSM.lookup("Minf");
  Info<< "Minf = " << Minf.value() <<  nl << endl;
  Ms = Minf;
}

//************************************************************//


Info<< "Reading thermodynamicProperties\n" << endl;

IOdictionary thermodynamicProperties
(
  IOobject
  (
    "thermodynamicProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
  )
);


dimensionedScalar R
(
  thermodynamicProperties.lookup("R")
);

dimensionedScalar Cv
(
  thermodynamicProperties.lookup("Cv")
);


dimensionedScalar Cp = R+Cv;

volScalarField gamma
(
  IOobject
  (
    "gamma",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("gamma", Cp/Cv)

);




volScalarField RS
(
  IOobject
  (
    "RS",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("RS", R)
);


volScalarField CvS
(
  IOobject
  (
    "CvS",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("CvS", Cv)
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
  IOobject
  (
    "transportProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
  )
);



dimensionedScalar mu0
(
  transportProperties.lookup("mu0")
);

volScalarField mu0S
(
  IOobject
  (
    "mu0S",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("mu0S", mu0)
);



//Dedner Diss. coeff
dimensionedScalar cd
(
  transportProperties.lookup("cd")
);






Info<< "Reading field p\n" << endl;
volScalarField p
(
  IOobject
  (
    "p",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);



Info<< "Reading field T\n" << endl;
volScalarField T
(
  IOobject
  (
    "T",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

volScalarField psi
(
  IOobject
  (
    "psi",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  1.0/(R*T)
);


volScalarField rho
(
  IOobject
  (
    "rho",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  psi*p
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
  IOobject
  (
    "U",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

volVectorField rhoU
(
  IOobject
  (
    "rhoU",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  rho*U
);


volScalarField rhoUx
(
  IOobject
  (
    "rhoUx",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  rhoU.component(vector::X)
);


volScalarField rhoUy
(
  IOobject
  (
    "rhoUy",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  rhoU.component(vector::Y)
);


volScalarField rhoUz
(
  IOobject
  (
    "rhoUz",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  rhoU.component(vector::Z)
);


Info<< "Reading field B components\n" << endl;


volScalarField Bx
(
  IOobject
  (
    "Bx",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);


volScalarField By
(
  IOobject
  (
    "By",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);



volScalarField Bz
(
  IOobject
  (
    "Bz",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);

//Assemble B vector
volVectorField B
(
  IOobject
  (
    "B",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  Bx*vector(1.0,0.0,0.0)
  + By*vector(0.0,1.0,0.0)
  + Bz*vector(0.0,0.0,1.0)
);


//*************************************************************
// Check if the case is magnetic or gas-dynamic only
// The flux method AUSMupPlus only works with GD

bool magnetic(false);
if (max(mag(B.internalField())) > 0.0)
{
  magnetic = true;
}
else if   (max(mag(B.boundaryField())) > 0.0)
{
  magnetic = true;
}

//*************************************************************
// Check if the case is 1D or multi-D
// Only works with grids aligned in x-direction)
// If the case is 1D Bx=const and no divB correction is needed

volVectorField OneD = fvc::reconstruct(mesh.magSf());

bool multiDim(true);
if (max(mag(OneD.internalField())) == 1.0)
{
  multiDim = false;
}
//*************************************************************


volScalarField pB
(
  IOobject
  (
    "pB",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);



volScalarField rhoE
(
  IOobject
  (
    "rhoE",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  p/(gamma-1.0) + rho*0.5*magSqr(U) + magSqr(B)*0.5/(mu0)
);



surfaceScalarField pos
(
  IOobject
  (
    "pos",
    runTime.timeName(),
    mesh
  ),
  mesh,
  dimensionedScalar("pos", dimless, 1.0)
);


surfaceScalarField neg
(
  IOobject
  (
    "neg",
    runTime.timeName(),
    mesh
  ),
  mesh,
  dimensionedScalar("neg", dimless, -1.0)
);




//Initialize the inviscid and magnetic fluxes
surfaceScalarField phi
(
  IOobject
  (
    "phi",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh.Sf() & linearInterpolate(rhoU)
);



surfaceScalarField Bf
(
  IOobject
  (
    "Bf",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh.Sf()/mesh.magSf() & linearInterpolate(B)
);

volScalarField divBf
(
  IOobject
  (
    "divBf",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  fvc::div(Bf)
);


surfaceScalarField Uf
(
  IOobject
  (
    "Uf",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh.Sf()/mesh.magSf() & linearInterpolate(U)
);



surfaceVectorField phiUpB
(
  IOobject
  (
    "phiUpB",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  phi*linearInterpolate(U)
  -Bf*linearInterpolate(B)/mu0*mesh.magSf()
  +linearInterpolate(p)*mesh.Sf()
);


surfaceVectorField phiUB
(
  IOobject
  (
    "phiUB",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  (Uf*linearInterpolate(B)
  -Bf*linearInterpolate(U))*mesh.magSf()
);


surfaceScalarField phiEpB
(
  IOobject
  (
    "phiEpB",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  Uf*linearInterpolate(rhoE + p + magSqr(B)/(2.0*mu0))*mesh.magSf()
);


//Co variables initialization
surfaceScalarField cfmax
(
  IOobject
  (
    "cfmax",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  Uf*mesh.magSf()
);
