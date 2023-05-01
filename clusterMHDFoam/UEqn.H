
//rhoUxEqn
solve(
      fvm::ddt(rhoUx)
      +fvc::div(phiUpB.component(vector::X))
      );

//rhoUyEqn
solve(
      fvm::ddt(rhoUy)
      +fvc::div(phiUpB.component(vector::Y))
      );

//rhoUzEqn
solve(
      fvm::ddt(rhoUz)
      +fvc::div(phiUpB.component(vector::Z))
      );

//Re-BUILD rhoU vector
rhoU = rhoUx*vector(1.0,0.0,0.0)
     + rhoUy*vector(0.0,1.0,0.0)
     + rhoUz*vector(0.0,0.0,1.0);

// Update internal field of U
U.dimensionedInternalField() =rhoU.dimensionedInternalField()
                             /rho.dimensionedInternalField();

// Correct boundary conditions of U
U.correctBoundaryConditions();
