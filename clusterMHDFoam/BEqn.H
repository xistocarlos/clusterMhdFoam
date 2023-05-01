
if(multiDim)
  {
    //BxEqn
    solve (
	   fvm::ddt(Bx)
	   +fvc::div(phiUB.component(vector::X))
	   );
  }
//ByEqn
solve(
      fvm::ddt(By)
      +fvc::div(phiUB.component(vector::Y))
      );

//BzEqn
solve(
      fvm::ddt(Bz)
      +fvc::div(phiUB.component(vector::Z))
      );


Bx.correctBoundaryConditions(); //Correct Bx BC
By.correctBoundaryConditions(); //Correct By BC
Bz.correctBoundaryConditions(); //Correct Bz BC

if(multiDim)
  {

    //Dedner hyperbolic/parabolic correction

    dimensionedScalar dT
      (
       "dT",
       dimensionSet(0,0,1,0,0,0,0),
       scalar(runTime.deltaT().value())
       );


    dimensionedScalar dx
      (
       "dx",
       dimensionSet(0,-1,0,0,0,0,0),
       scalar(max(mesh.surfaceInterpolation::deltaCoeffs()).value())
       );

    dimensionedScalar ch = CoNum/(dT*dx);
    dimensionedScalar cp = Foam::sqrt(-dT*sqr(ch)/Foam::log(cd));

    solve
      (
       fvm::ddt(pB)
       +sqr(ch)*fvc::div(Bf*mesh.magSf())
       +fvm::Sp(sqr(ch)/sqr(cp),pB)
       );

    pB.correctBoundaryConditions();
  }

//Rebuild B vector
B.dimensionedInternalField() = Bx.dimensionedInternalField()*vector(1.0,0.0,0.0)
  + By.dimensionedInternalField()*vector(0.0,1.0,0.0)
  + Bz.dimensionedInternalField()*vector(0.0,0.0,1.0);

B.boundaryField() = Bx.boundaryField()*vector(1.0,0.0,0.0)
  + By.boundaryField()*vector(0.0,1.0,0.0)
  + Bz.boundaryField()*vector(0.0,0.0,1.0);
