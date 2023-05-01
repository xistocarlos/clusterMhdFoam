/*
=============================================================================
AUSMpw-MHD FLUX START

This subroutine computes the AUSMpw-MHD flux for the MHD equations
in the direction, S=[Sx,Sy,Sz].
Conservative form of the MHD equations:

     dU/dt + dF/dx + dG/dy + dH/dz = 0


     Fn = F*Sx + G*Sy + H*Sz = | rho*Uf                  |
                               | rho*Uf*u - BxBf + pG*Sx |
                               | rho*Uf*v - ByBf + pG*Sy |
                               | rho*Uf*w - BzBf + pG*Sz |
                               |    Uf*Bx - Bf*u         |
                               |    Uf*By - Bf*v         |
                               |    Uf*Bz - Bf*w         |
                               | rho*Uf*H - Bf(U & B)    |  

(Uf =  u*Sx +  v*Sy +  w*Sz)
(Bf = Bx*Sx + By*Sy + Bz*Sz)

The flux function is given by:

F(i+1/2) = cf*(M4(+)*PhiL+M4(-)*PhiR) + (P5(+)*PL+P5(-)PR) + 0.5*(PhiBL+PhiBR)


     Phi = | rho   |      P = |    0      |     PhiB = |    0   |    
           | rho*u |          |  Sx*pG    |            | -Bx*Bf |
           | rho*v |          |  Sy*pG    |            | -By*Bf |
           | rho*w |          |  Sz*pG    |            | -Bz*Bf |
           |  Bx   |          | -Bf*Ux    |            |    0   |             
           |  By   |          | -Bf*Uy    |            |    0   | 
           |  Bz   |          | -Bf*Uz    |            |    0   | 
           |  Ht   |          | -Bf*(U&B) |            |    0   | 



The M4(+)(-) are modified by the introuction of global pressure-based weigth 
functions.


Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
         flux description. Therefore do not stress to much if you find a pB*Sf below. 



References   
(1) Han, Sang-Hoon, Jeong-Il Lee, and Kyu Hong Kim. "Accurate and robust pressure
    weight advection upstream splitting method for magnetohydrodynamics equations"
    AIAA journal 47.4 (2009): 970-981.
(2) Xisto, C. M., Páscoa, J. C. , Oliveira, P .J. (2013), “A pressure-based method
    with AUSM-type fluxes for MHD flows at arbitrary Mach numbers”, International 
    Journal for Numerical Methods in Fluids, Vol. 72:11, pages 1165-1182, 
    doi:10.1002/fld.3781.
(3) Xisto, C. M., Páscoa, J. C., Oliveira, P. J. (2014), “A pressure-based high
    resolution numerical method for resistive MHD”, Journal of Computational
    Physics, Vol. 275, pages 323-345, doi:10.1016/j.jcp.2014.07.009.


Authors: Carlos Xisto (xisto@ubi.pt)
Date: ..-..-2014
=============================================================================
*/


void AUSMpwFlux (
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
		 scalar &Bf, scalar &rhoflux,  
		 vector &rhoUflux, vector &rhoBflux, 
		 scalar &rhoEflux, scalar &cfmax
		 )
 
{

  const scalar beta = 1.0/8.0;			                       
  const scalar alpha = 3.0/16.0;
  const scalar mu0 = 1.0; //temos de mudar


  //================== Direction Spliting Uf = U & Sf | Bf = B & Sf =======================
  const vector Sf_unit = Sf /magSf;

  const scalar BLf = BL  & Sf_unit;       
  const scalar BRf = BR  & Sf_unit;     
           
  const scalar ULf = UL & Sf_unit;       
  const scalar URf = UR & Sf_unit;     
                    											
  //================= CALCULATE magSqr(B) and magSqr(U) ==================================== 
  const scalar BLsqr = magSqr(BL); 
  const scalar BRsqr = magSqr(BR); 

 
  //======================================================================================== 
  //    CALCULATE NUMERICAL SOUND SPEED (Eq(29, liou2006)) and FAST MAGNETOSONIC SPEED
  //======================================================================================== 

  //------------ TOTAL ENTALPHY ------------------------------------------------------------
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




  //========================================================================================== 
  //			CALCULATE P(5) and M(4) FUNCTIONS
  //========================================================================================== 

 
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

  /*
    const scalar M4pos = ((magML  >= 1.0) ? M1pos   
    : (0.25*sqr(ML+1.0)+beta*sqr(sqr(ML)-1.0)));
                            
    const scalar M4neg = ((magMR  >= 1.0) ? M1neg   
    : (-0.25*sqr(MR-1.0)-beta*sqr(sqr(MR)-1.0)));

    const scalar P5pos = ((magML  >= 1.0) ? M1pos/(ML+SMALL)
    : (0.25*sqr(ML+1.0)*(2.0-ML)+alpha*ML*sqr(sqr(ML)-1.0)));

    const scalar P5neg = ((magMR  >= 1.0) ? M1neg/(MR+SMALL)
    : (0.25*sqr(MR-1.0)*(2.0-MR)-alpha*MR*sqr(sqr(MR)-1.0)));
  */

  /*
    ========================================================================================== 
    AUSMPW-HAN Weigth-functions FOR MHD   
    ========================================================================================== 
  */


  const scalar pLt = pL + BLsqr/(2.0*mu0);  
  const scalar pRt = pR + BRsqr/(2.0*mu0);

  const scalar w = 1.0-Foam::pow(min(pLt/pRt, pRt/pLt),3.0); //Eqn(15) Han 2009

  const scalar ps = (P5pos*pLt + P5neg*pRt);  //Eqn(13) Han 2009 

  const scalar fR = ((mag(ps)  != 0.0) ? ((pRt)/(ps)-1.0) : 0.0); //Eqn(14) Han 2009 
  const scalar fL = ((mag(ps)  != 0.0) ? ((pLt)/(ps)-1.0) : 0.0);
			

	
  const scalar Mn=M4pos + M4neg;
	
  const scalar M4posN = ((Mn  >= 0.0) ? (M4pos + M4neg*((1.0-w)*(1.0+fR)-fL)) 
			 : (M4pos*w*(1.0+fL)));


                             
  const scalar M4negN = ((Mn  >= 0.0) ? (M4neg*w*(1.0+fR))
			 : (M4neg + M4pos*((1.0-w)*(1.0+fL)-fR)));


  //=====================================================================================




  //Magnetic Flux
  Bf = (BLf + BRf)/2.0;      

  //CONTINUITY FLUX 
  rhoflux = cf*(M4posN*rhoL + M4negN*rhoR)*magSf;


  //MOMENTUM FLUX  

  const scalar pGf = (pL+BLsqr/(2.0*mu0))*P5pos + (pR+BRsqr/(2.0*mu0))*P5neg;
  const vector BfB = (Bf*BL +  Bf*BR)/(2.0*mu0);


  rhoUflux = cf*(M4posN*rhoUL + M4negN*rhoUR)*magSf
    - BfB*magSf + pGf*Sf; 
             
  //INDUCTION FLUX  
  const scalar  pBf = 0.5*(pBL + pBR);
  const vector  UfB = cf*(M4posN*BL+M4negN*BR);
  const vector  BfU = (P5pos*( Bf*UL)+P5neg*( Bf*UR));
  rhoBflux = (UfB-BfU)*magSf + pBf*Sf;


  //ENERGY FLUX 
  const scalar HtplusBL = HtL + BLsqr/(2.0*mu0);
  const scalar HtplusBR = HtR + BRsqr/(2.0*mu0);
  
  const scalar BUL = UL & BL;
  const scalar BUR = UR & BR; 
  const scalar BfBU = P5pos*( Bf*BUL)/mu0 + P5neg*( Bf*BUR)/mu0;
   
  rhoEflux = (cf*(M4posN*HtplusBL+M4negN*HtplusBR) - BfBU)*magSf;

  //Co variables
  const scalar Uf = cf*(M4posN+M4negN);

  cfmax = (mag(Uf)+cf)*magSf;
  
  //============================================================================
  //                   AUSMpw-MHD END                       
  //============================================================================

}
