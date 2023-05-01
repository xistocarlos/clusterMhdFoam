
if(multiDim)
  {
    //BxEqn
    Bx=B0.component(vector::X)
      - RK4alpha[i]*runTime.deltaT()
      *fvc::div(phiUB.component(vector::X));
  }
//ByEqn	 
By=B0.component(vector::Y)
  - RK4alpha[i]*runTime.deltaT()
  *fvc::div(phiUB.component(vector::Y));

//BzEqn	 
Bz=B0.component(vector::Z)
  - RK4alpha[i]*runTime.deltaT()
  *fvc::div(phiUB.component(vector::Z));


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
    

    pB=pB0-RK4alpha[i]*runTime.deltaT()*
      (
       sqr(ch)*fvc::div(Bf*mesh.magSf())
       +(sqr(ch)/sqr(cp))*pB
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



