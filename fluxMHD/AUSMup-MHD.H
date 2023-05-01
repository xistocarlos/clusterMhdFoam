/*
  =============================================================================
  AUSM+-up FLUX START

  This subroutine computes the AUSM+-up-MHD flux for the MHD equations
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


  The mass flux function is given by:
                   |rhoL if Mf>0
  F(i+1/2) = cf*Mf<               
                   |rhoR otherwise

  where cf = min(cL,cR) (Now is computed acording to the fast magnetosonic speed) 
  and 
  Mf = M4pos + M4neg + Mp.
  M4pos/neg are mach numeber interpolation functions and
  Mp is a pressure diffusion term introduced to enhace calculations at low Mach.

  The pressure flux function is given by:

  pf = P5pos*pL + P5neg*pR + Pu;

  Where P5pos/neg are again Mach number interpolation functions and
  Pu is a velocity diffusion terms 

  Modification for MHD:
  Replacment of acoutic wave: 

  |Uf|+af

  by the magnetosonic fast wave:

  |Uf|+cf

  Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
  flux description. Therefore do not stress to much if you find a pB*Sf below. 

  References
  (1) Meng-Sing Liou, "A sequel to AUSM, Part II: AUSM+-up for all speeds"
  JCP, 2006.

  (2) Xisto, C. M., Páscoa, J. C. , Oliveira, P. J., Nicolini, D. A. (2012), 
  An hybrid pressure­ density based algorithm for the Euler equations at 
  all Mach number regimes, International Journal for Numerical Methods 
  in Fluids Vol. 70:8, pages 961-976, doi:10.1002/fld.2722.



  Authors: Carlos Xisto (xisto@ubi.pt)
  Date: 13-02-2015
  =============================================================================
*/

void AUSMupMHDFlux (
		    const vector &Sf,
		    const scalar &magSf,
		    const scalar &rhoL,  const scalar &rhoR, 
		    const vector &rhoUL, const vector &rhoUR, 
		    const vector &UL,    const vector &UR, 
		    const scalar &pL,    const scalar &pR, 
		    const scalar &rhoEL, const scalar &rhoER, 
		    const scalar &gammaL, const scalar &gammaR, 
		    const scalar &Ms,
		    const scalar &pBL,   const scalar &pBR, 
		    const vector &BL,    const vector &BR, 
		    scalar &Bf, scalar &rhoflux,  
		    vector &rhoUflux, vector &rhoBflux, 
		    scalar &rhoEflux, scalar &cfmax
		    )
{


  const scalar Ku = 0.75; 
  const scalar Kp = 0.25;
  const scalar beta = 1.0/8.0; 
  const scalar mu0 = 1.0; //temos de mudar


  //================== Direction Spliting Uf = U & Sf    =======================
  const vector Sf_unit = Sf /magSf;
      
  const scalar    ULf = UL & Sf_unit;       
  const scalar    URf = UR & Sf_unit;  

  const scalar BLf = BL  & Sf_unit;       
  const scalar BRf = BR  & Sf_unit;   

  //================= CALCULATE magSqr(B) and magSqr(U) ==================================== 
  const scalar BLsqr = magSqr(BL); 
  const scalar BRsqr = magSqr(BR);   
 	
  //============================================================================
  //    CALCULATE NUMERICAL SOUND SPEED (Eq(29, liou2006)) 
  //=============================================================================

  //------------ TOTAL ENTALPHY -------------------------------------------------

  const scalar HtL = rhoEL + pL; 
  const scalar HtR = rhoER + pR;
 
  const scalar Htf = min(HtL/rhoL, HtR/rhoR); 

  const scalar gammaF = 0.5*(gammaL+gammaR); 
  
	
  //------------NUMERICAL SOUND SPEED---------------------------------------------------------
  const scalar acritF = Foam::sqrt(mag((2.0*Htf*(gammaF-1.0))/(gammaF+1.0))); 

  const scalar    aL = sqr(acritF)/(max(acritF, ULf));   //Eq(30, liou2006)
  const scalar    aR = sqr(acritF)/(max(acritF,-URf));   

 

  //-------FAST MAGNETOSONIC SPEED------------------------------------------------------------
  const scalar bfL = BLsqr/(rhoL*mu0); //alfven speed
  const scalar bfR = BRsqr/(rhoR*mu0); //alfven speed

  const scalar caL = BLf/Foam::sqrt(rhoL*mu0); //normal alfven speed
  const scalar caR = BRf/Foam::sqrt(rhoR*mu0); 

  const scalar cfL = Foam::sqrt((1.0/2.0) * (sqr(aL) + bfL + Foam::sqrt(sqr(sqr(aL)+ bfL) 
									- 4.0*sqr(aL)*sqr(caL) )));//left value for magnetosonic fast speed
  const scalar cfR = Foam::sqrt((1.0/2.0) * (sqr(aR) + bfR + Foam::sqrt(sqr(sqr(aR)+ bfR) 
									- 4.0*sqr(aR)*sqr(caR) )));//rigth value for magnetosonic fast speed

  const scalar   cf = Foam::min(cfL, cfR);

  //--------------------CALCULATE MACH NUMBER Eq(16, liou2006)--------------------------------
  const scalar  ML = ULf/cf;
  const scalar  MR = URf/cf; 



  //Calcular precondicionamento Eq(11, liou2006)
  const scalar Mmed = Foam::sqrt((sqr(ML)+sqr(MR))*0.5);

  //eq(14 e 71, Liou2006)
  const scalar Mref = Foam::sqrt(min(1.0, max(sqr(Mmed), sqr(Ms)))); 
  const scalar fpre = Mref*(2.0-Mref); //eq(14, Liou2006) 

  //=============================================================================
  //			CALCULATE P(5) and M(4) FUNCTIONS
  //=============================================================================
                 
  const scalar alpha = (3.0/16.0)*(-4.0+5*sqr(fpre));             

  const scalar   M1pos = (ML+mag(ML))/2.0; 
  const scalar   M1neg = (MR-mag(MR))/2.0; 


  const scalar magMR = mag(MR);
  const scalar magML = mag(ML);

  const scalar M4pos = Foam::pos(magML - 1.0)*M1pos 
    + Foam::neg(magML - 1.0)*(
			      0.25*sqr(ML+1.0)+beta*sqr(sqr(ML)-1.0)
			      );

  const scalar M4neg = Foam::pos(magMR - 1.0)*M1neg 
    + Foam::neg(magMR - 1.0)*( 
			      -0.25*sqr(MR-1.0)-beta*sqr(sqr(MR)-1.0)
			       );
                            
  const scalar P5pos = Foam::pos(magML - 1.0)*(M1pos/(ML+SMALL)) 
    + Foam::neg(magML - 1.0)*(
			      0.25*sqr(ML+1.0)*(2.0-ML)+alpha*ML*sqr(sqr(ML)-1.0)
			      );
                            
  const scalar P5neg = Foam::pos(magMR - 1.0)*(M1neg/(MR+SMALL)) 
    + Foam::neg(magMR - 1.0)*(
			      0.25*sqr(MR-1.0)*(2.0+MR)-alpha*MR*sqr(sqr(MR)-1.0)
			      );


  //=============================================================================
  //                            Calcular os Termos dissipativos
  //=============================================================================
  const scalar rhoM = (rhoL + rhoR)/2;

  const scalar   Mp = -(Kp/fpre)*max((1.0-sqr(Mmed)),0.0) 
    * ((pR - pL)/(rhoM*sqr(cf))); //(Eq73)

  const scalar   Pu = -Ku*P5pos*P5neg*(rhoL+rhoR)
    *(fpre*(cf))*((URf-ULf));   //termo difusivo (Eq75) 


  //=============================================================================
  //                                  FLUX COMPUTATION
  //=============================================================================

  const scalar Uf = cf*(M4pos+M4neg+Mp);      



  //Magnetic Flux
  Bf = (BLf + BRf)/2.0; 

  //CONTINUITY FLUX 
  rhoflux =  0.5*(Uf*(rhoL+rhoR)-mag(Uf)*(rhoR-rhoL))*magSf;


  //MOMENTUM FLUX  

  const scalar pGf = (pL+BLsqr/(2.0*mu0))*P5pos + (pR+BRsqr/(2.0*mu0))*P5neg + Pu;
  const vector BfB = (Bf*BL +  Bf*BR)/(2.0*mu0);

  rhoUflux = 0.5*(Uf*(rhoUL+rhoUR)-mag(Uf)*(rhoUR-rhoUL))*magSf
    - BfB*magSf + pGf*Sf;
 
  //INDUCTION FLUX  
  const scalar  pBf = 0.5*(pBL + pBR);
  const vector  BfU = (P5pos*( Bf*UL)+P5neg*( Bf*UR));

  rhoBflux = 0.5*(Uf*(BL+BR)-mag(Uf)*(BR-BL))*magSf 
    - BfU*magSf + pBf*Sf;


  //ENERGY FLUX 
  const scalar HtplusBL = HtL + BLsqr/(2.0*mu0);
  const scalar HtplusBR = HtR + BRsqr/(2.0*mu0);
  
  const scalar BUL = UL & BL;
  const scalar BUR = UR & BR; 
  const scalar BfBU = P5pos*( Bf*BUL)/mu0 + P5neg*( Bf*BUR)/mu0;
   

  rhoEflux = 0.5*(Uf*(HtplusBL+HtplusBR)-mag(Uf)*(HtplusBR-HtplusBL))*magSf
    - BfBU*magSf; 

  //Co variables
  cfmax = (mag(Uf)+cf)*magSf;

}                                          
//============================================================================= 
//                                   AUSM+ FLUX END
//=============================================================================

