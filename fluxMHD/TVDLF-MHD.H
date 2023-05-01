/*
  =============================================================================
  TVDLF FLUX START

  This subroutine computes the TVDLF flux for the MHD equations
  in the direction, S=[Sx,Sy,Sz].
  Conservative form of the MHD equations:

  dU/dt + dF/dx + dG/dy + dH/dz = 0


  Fn = F*Sx + G*Sy + H*Sz = | rho*Uf                  |
                            | rho*Uf*u - BxBf + pG*Sx |
                            | rho*Uf*v - ByBf + pG*Sy |
                            | rho*Uf*w - BzBf + pG*Sz |
                            | Uf*Bx - Bf*u + Psi*Sx   |
                            | Uf*By - Bf*v + Psi*Sy   |
                            | Uf*Bz - Bf*w + Psi*Sz   |
                            | rho*Uf*H - Bf(U & B)    |  

  (Uf =  u*Sx +  v*Sy +  w*Sz)
  (Bf = Bx*Sx + By*Sy + Bz*Sz)

  The flux function is given by: 

  F(i+1/2) = 1/2*[F(UL) + F(UR) cmax*(UR-UL)]

  where: cmax is the local maximumal physical propagation speed (fast magnetosonic):
    
  cmax= max(|ULf|+cL,|URf|+cR) or cmax=1/2*[mag(ULf+URf)+(cL+cR)]; 

  This method is diffusive, but extremely robust. It can be emplyed for 
  computing zero-beta plasma flow (in MHD version).

  Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
  flux description. Therefore do not stress to much if you find a pB*Sf below. 

  References:
  (1) Toth and Odstrcil, (1996)"Comparison of Some Flux Correced Transport
  and Total Variation Diminishing Numerical Schemes for Hydrodynamics and 
  Magnetohydrodynamics Problems", Journal of Computational Physics 128,82

  Author: Carlos Xisto
  date: 1-1-2015
  =============================================================================
*/

void TVDLFMHDFlux (
		   const vector &Sf,
		   const scalar &magSf,
		   const scalar &rhoL,  const scalar &rhoR, 
		   const vector &rhoUL, const vector &rhoUR, 
		   const vector &UL,    const vector &UR, 
		   const scalar &pL,    const scalar &pR, 
		   const scalar &rhoEL, const scalar &rhoER, 
		   const scalar &gammaL, const scalar &gammaR, 
		   const scalar &pBL,   const scalar &pBR, 
		   const vector &BL,    const vector &BR, 
		   scalar &Bf, scalar &rhoflux, vector &rhoUflux, 
		   vector &rhoBflux,   scalar &rhoEflux, scalar &cfmax
		   )
{


  const scalar mu0 = 1.0; //tem de se mudar


  //================== Direction Spliting Uf = U & Sf | Bf = B & Sf =======================
  const vector Sf_unit = Sf/magSf;

  const scalar BLf = BL  & Sf_unit;       
  const scalar BRf = BR  & Sf_unit;     

  const scalar ULf = UL & Sf_unit;       
  const scalar URf = UR & Sf_unit;  

 

  //-------FAST MAGNETOSONIC SPEED------------------------------------------------------------


  const scalar BLsqr = magSqr(BL); 
  const scalar BRsqr = magSqr(BR);

  const scalar    aL = Foam::sqrt(gammaL*pL/rhoL);   
  const scalar    aR = Foam::sqrt(gammaR*pR/rhoR); 
 
  const scalar bfL = BLsqr/(rhoL*mu0); //alfven speed
  const scalar bfR = BRsqr/(rhoR*mu0); //alfven speed

  const scalar caL = BLf/Foam::sqrt(rhoL*mu0); //normal alfven speed
  const scalar caR = BRf/Foam::sqrt(rhoR*mu0); 


  const scalar cfL = Foam::sqrt((1.0/2.0) * (sqr(aL) + bfL + Foam::sqrt(sqr(sqr(aL)+ bfL) 
									- 4.0*sqr(aL)*sqr(caL) )));//left value for magnetosonic fast speed
  const scalar cfR = Foam::sqrt((1.0/2.0) * (sqr(aR) + bfR + Foam::sqrt(sqr(sqr(aR)+ bfR) 
									- 4.0*sqr(aR)*sqr(caR) )));//rigth value for magnetosonic fast speed


  //const scalar    lambda1 = mag(ULf)+cfL;   
  //const scalar    lambda2 = mag(URf)+cfR;   
  //const scalar    lambdaMax = max(lambda1,lambda2);

  const scalar    lambdaMax = 0.5*mag(ULf+URf)+0.5*(cfL+cfR);

  //Magnetic Flux
  Bf = (BLf + BRf)/2.0; 
  /*
    It is necessary for the pB equation      
  */

  //CONTINUITY FLUX 
  rhoflux = 0.5*(URf*rhoL  + ULf*rhoR  - lambdaMax*(rhoR-rhoL))*magSf;


  //MOMENTUM FLUX  

  const scalar pGf = pL+BLsqr/(2.0*mu0) + pR+BRsqr/(2.0*mu0);
  const vector BfB = (BLf*BL + BRf*BR)/(mu0);


  rhoUflux = 0.5*( ULf*rhoUL + URf*rhoUR  - BfB + pGf*Sf_unit 
		   - lambdaMax*(rhoUR-rhoUL))*magSf;
  
  
  //INDUCTION FLUX
  const scalar pBf = (pBL + pBR);
  const vector BfU = (BLf*UL)+(BRf*UR);
  rhoBflux = 0.5*(ULf*BL + URf*BR - BfU 
		  + pBf*Sf_unit - lambdaMax*(BR-BL))*magSf;  
  



  //ENERGY FLUX 
  const scalar BUL = UL & BL;
  const scalar BUR = UR & BR; 
  const scalar BfBU = (BLf*BUL)/mu0 + (BRf*BUR)/mu0; 
  
  const scalar pGL = pL + BLsqr/(2.0*mu0);
  const scalar pGR = pR + BRsqr/(2.0*mu0);
   
  rhoEflux = 0.5*(ULf*(rhoEL+pGL) + URf*(rhoER+pGR) 
		  - BfBU - lambdaMax*(rhoER-rhoEL))*magSf;

  
  
  cfmax =lambdaMax*magSf;
  
  //============================================================================
  //                   TVDLF END                       
  //============================================================================

}
