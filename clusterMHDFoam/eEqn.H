    
solve
(
 fvm::ddt(rhoE)
 + fvc::div(phiEpB)       
 );


T.dimensionedInternalField() = 
  (
   rhoE.dimensionedInternalField()
   -0.5*rho.dimensionedInternalField()
   *magSqr(U.dimensionedInternalField())
   -0.5*magSqr(B.dimensionedInternalField())
   /(mu0S.dimensionedInternalField())
   )
  /CvS.dimensionedInternalField()
  /rho.dimensionedInternalField();

T.correctBoundaryConditions();  


//Update psi  field
psi.dimensionedInternalField()= 
  1.0/(T.dimensionedInternalField()
       *RS.dimensionedInternalField());   

psi.boundaryField()= 
  1.0/(T.boundaryField()
       *RS.boundaryField());



// Calculate internal pressure field    
p.dimensionedInternalField() =
  rho.dimensionedInternalField()
  /psi.dimensionedInternalField();

// Update boundary conditions for p          
p.correctBoundaryConditions(); 

//Update boundary field of rho           
rho.boundaryField() = psi.boundaryField()*p.boundaryField();



//Update boundary field of rhoU 
rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
rhoUx.boundaryField() =  rhoU.boundaryField().component(vector::X);
rhoUy.boundaryField() =  rhoU.boundaryField().component(vector::Y);
rhoUz.boundaryField() =  rhoU.boundaryField().component(vector::Z);



// Update boundary field of rhoE
rhoE.boundaryField() =
  rho.boundaryField()*
  (
   T.boundaryField()*CvS.boundaryField() 
   + 0.5*magSqr(U.boundaryField())
   )
  + 0.5*magSqr(B.boundaryField())/(mu0S.boundaryField());







