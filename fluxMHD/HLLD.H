/*
  =============================================================================
  HLLD FLUX START
  This subroutine computes the HLLD flux for the MHD equations
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


           |Fn(UL)  if  SL>0,  
           |F(UL*)  if  SL<0<SL*,
           |F(UL**) if SL*<0<SM,
  FHLLD = <
           |F(UR**) if  SM<0<SR*
           |F(UR*)  if SR*<0<SR,
           |Fn(UR)  if  Sr<0,


  SL  = min(Uf-af, UfTile-aTilde);
  SL* = 
  SR* =
  SR  = max(Uf+af, UfTile+aTilde);

  Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
  flux description. Therefore do not stress to much if you find a pB*Sf below. 

  
  Refs:
  (1) P. Batten, M.A. Leschziner, U.C. Goldberg, Average-State Jacobians and
  Implicit Methods for Compressible Viscous and Turbulent Flows, Journal of 
  Computational Physics, Volume 137, Issue 1, October 1997, Pages 38-78, 
  ISSN 0021-9991, http://dx.doi.org/10.1006/jcph.1997.5793

  (2) Miyoshi and Kusano (2005) "A multi-state HLL aproximate Riemann solver for
  ideal MHD", JCP, vol. 208, pp. 315-344.


  Authors: Carlos Xisto (xisto@ubi.pt)
  Date: work in progress
  =============================================================================
*/


void HLLDFlux (
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

  const scalar mu0 = 1.0; //temos de mudar

  //--------------------- Direction Spliting Uf = U & Sf | Bf = B & Sf ---------
  const vector Sf_unit = Sf /magSf;
           
  const scalar ULf = UL & Sf_unit;       
  const scalar URf = UR & Sf_unit;     

  const scalar BLf = BL & Sf_unit;       
  const scalar BRf = BR & Sf_unit;   
                    											
  //-------FAST MAGNETOSONIC SPEED----------------------------------------------


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

  //-------Signal speeds--------------------------------------------------------

  const scalar SL = min(ULf,URf) - max(cfL,cfR); //Eq. 67 
  const scalar SR = max(ULf,URf) + max(cfL,cfR);

  //const scalar SL = min(ULf-cfL, URf-cfR); //Eq. 12 
  //const scalar SR = max(ULf+cfL, URf+cfR);

  const scalar pGL = pL+BLsqr/(2.0*mu0);
  const scalar pGR = pR+BRsqr/(2.0*mu0);  

  const scalar SM = (rhoL*ULf*(ULf-SL) - rhoR*URf*(URf-SR) + pGL - pGR)/
    (rhoL*(ULf-SL)-rhoR*(URf-SR));//Eq.38

  const scalar SLstar = SM - caL;
  const scalar SRstar = SM + caR;


  //------Compute pressure in star region from the right and left sides-------------------
  const scalar pGstarR = pGR + rhoR*(URf-SR)*(URf-SM);//Eq.23
  const scalar pGstarL = pGL + rhoL*(ULf-SL)*(ULf-SM);//Eq.23

  //use pGStarR for pGStar, as in theory, pGL*=PGL**=pGR**=pGR*=pG*
  const scalar pGstar = pGstarR;

  //additionally one can specify
  //const scalar pGstar = ((SR-URf)*rhoR*pGL-(SL-ULf)*rhoL*pGR+rhoL*rhoR*(SR-URf)
  //                      *(SL-ULf)*(URf-ULf))/((SR-URf)*rhoR-(SL-ULf)*rhoL); //Eq.41

  // give a warning if this is not the case
  if ( mag(pGstarR-pGstarL) > 1e-6 )
    {
      Info << "WARNING, mag(pGStarRight-pGStarLeft) > SMALL " << endl;
    }

  /*
    The following implementation is similar to variation HLLC-L (only two-states)
    With the conistent formulation of Mignone (see Jets from Young Stars pp.84)
  */

  const scalar rhoLstar = rhoL*(SL-ULf)/(SL-SM);
  const scalar rhoRstar = rhoR*(SR-URf)/(SR-SM);

  //const scalar BfHll=(SR*BRf-SL*BLf)/(SR-SL); //eq(33) em Li(2005) 

  //const vector ULstar = UL-(Sf_unit*(pGL-pGstar)*(SL-SM)+BLf*BL*(SM-ULf))
  //                        /(rhoL*(SL-ULf)*(SL-SM)-sqr(BLf));//Eq. 46 in Miyoshi for 3D
  //const vector URstar = UR-(Sf_unit*(pGR-pGstar)*(SR-SM)+BRf*BR*(SM-URf))
  //                        /(rhoR*(SR-URf)*(SR-SM)-sqr(BRf));//Eq. 46 in Miyoshi for 3D

  //const vector BLstar = (BL*(rhoL*sqr(SL-ULf)-sqr(BLf))+BLf*Sf_unit*(pGL-pGstar))
  //                        /(rhoL*(SL-ULf)*(SL-SM)-sqr(BLf));//Eq. 47 in Miyoshi for 3D
  //const vector BRstar = (BR*(rhoR*sqr(SR-URf)-sqr(BRf))+BRf*Sf_unit*(pGR-pGstar))
  //                        /(rhoR*(SR-URf)*(SR-SM)-sqr(BRf));//Eq. 47 in Miyoshi for 3D



  const vector Bhll = (SR*BR-SL*BL + (ULf*BL - BLf*UL) 
		       - (URf*BR - BRf*UR))/(SR-SL);

  const vector rhoULstar = 1.0/(SL-SM)*(
					rhoUL*(SL-ULf)
					-(BLf*Bhll-BLf*BL)
					-Sf_unit*(pGL-pGstar)
					);

  const vector rhoURstar = 1.0/(SR-SM)*(
					rhoUR*(SR-URf)
					-(BRf*Bhll-BRf*BR)
					-Sf_unit*(pGR-pGstar)
					);

  const vector FBhll = 1.0/(SR-SL)*(
				    SR*(ULf*BL - BLf*UL)
				    -SL*(URf*BR - BRf*UR)
				    +SR*SL*(BR-BL)
				    );

  const vector Uhll = (Bhll*SM-FBhll)/stabilise(Bf,SMALL); //Eq. 43 do Mignone (J. f. y. stars)



  //const vector Uhll = (SR*UR-SL*UL + (ULf*rhoUL + pGL*Sf_unit - BLf*BL/mu0)
  //                                 - (URf*rhoUR + pGR*Sf_unit - BRf*BR/mu0))/(SR-SL);

  const scalar rhoELstar = (rhoEL*(SL-ULf)-pGL*ULf+pGstar*SM+BLf*((UL&BL) - (Uhll&Bhll)))
    /(SL-SM);
  const scalar rhoERstar = (rhoER*(SR-URf)-pGR*URf+pGstar*SM+BRf*((UR&BR) - (Uhll&Bhll)))
    /(SR-SM);



  if ( pos(SL) )
    {
      // compute FL
      rhoflux  = (ULf*rhoL)*magSf;
      rhoUflux = (ULf*rhoUL + pGL*Sf_unit - BLf*BL/mu0)*magSf;
      rhoBflux = (ULf*BL - BLf*UL)*magSf;
      rhoEflux = (ULf*(rhoEL+pGL) - BLf*((UL & BL))/mu0)*magSf;
    }
  else if ( pos(SM) )
    {
      // compute FL*
      rhoflux  = (SM*rhoLstar)*magSf;
      rhoUflux = (rhoULstar*SM + pGstar*Sf_unit - BLf*Bhll/mu0)*magSf;
      rhoBflux = (Bhll*SM - BLf*Uhll)*magSf;            
      rhoEflux = ((rhoELstar + pGstar)*SM - BLf*((Uhll & Bhll))/mu0)*magSf;
    }
  else if ( pos(SR) )
    {
      // compute FR*
      rhoflux  = (SM*rhoRstar)*magSf;
      rhoUflux = (rhoURstar*SM + pGstar*Sf_unit - BRf*Bhll/mu0)*magSf;
      rhoBflux = (Bhll*SM - BRf*Uhll)*magSf;            
      rhoEflux = ((rhoERstar + pGstar)*SM - BRf*((Uhll & Bhll))/mu0)*magSf;
    }
  else if ( neg(SR) )
    {
      // compute FR
      rhoflux  = (URf*rhoR)*magSf;
      rhoUflux = (URf*rhoUR + pGR*Sf_unit - BRf*BR/mu0)*magSf;
      rhoBflux = (URf*BR - BRf*UR)*magSf;
      rhoEflux = (URf*(rhoER+pGR) - BRf*((UR & BR))/mu0)*magSf;
    }
  else
    {
      Info << "Error in HLLD Riemann solver" << endl;
    }




 

  cfmax = max(mag(ULf)+cfL, mag(URf)+cfR)*magSf;

}


