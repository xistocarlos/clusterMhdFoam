
//rhoUxEqn
rhoUx=rhoU0.component(vector::X)
    -RK4alpha[i]*runTime.deltaT()
    *fvc::div(phiUpB.component(vector::X));


//rhoUyEqn
rhoUy=rhoU0.component(vector::Y)
     -RK4alpha[i]*runTime.deltaT()
     *fvc::div(phiUpB.component(vector::Y));

//rhoUzEqn
rhoUz=rhoU0.component(vector::Z)
     -RK4alpha[i]*runTime.deltaT()
     *fvc::div(phiUpB.component(vector::Z));

//Re-BUILD rhoU vector
rhoU = rhoUx*vector(1.0,0.0,0.0)
     + rhoUy*vector(0.0,1.0,0.0)
     + rhoUz*vector(0.0,0.0,1.0);

// Update internal field of U
U.dimensionedInternalField() =rhoU.dimensionedInternalField()
                             /rho.dimensionedInternalField();

// Correct boundary conditions of U
U.correctBoundaryConditions();
