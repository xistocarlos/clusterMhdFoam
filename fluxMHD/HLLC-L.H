/*
  =============================================================================
  HLLC-L FLUX START
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
  FHLLC-L = <
             |F(UR*)  if SR*<0<SR,
             |Fn(UR)  if  Sr<0,


  SL  = min(Uf-af, UfTile-aTilde);
  SR  = max(Uf+af, UfTile+aTilde);

  Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
  flux description. Therefore do not stress to much if you find a pB*Sf below. 

  
  Refs:
  (1) Shengtai Li, An HLLC Riemann solver for magneto-hydrodynamics, 
  Journal of Computational Physics, Volume 203, Issue 1, 10 February 2005,
  Pages 344-357, ISSN 0021-9991, http://dx.doi.org/10.1016/j.jcp.2004.08.020.

  (2) Massaglia, S. et al. (Eds.), Jets From Young Stars III, Lect. 
  Notes Phys. 754 (Springer, Berlin Heidelberg 2008),
  DOI 10.1007/ 978-3-540-76967-5

  TODO: Problems in defining UHll

  Authors: Carlos Xisto (xisto@ubi.pt)
  Date: 15/02/2015
  =============================================================================
*/


void HLLCLFlux (
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

  //const scalar SL = min(ULf,URf) - max(cfL,cfR); //Eq. 67 
  //const scalar SR = max(ULf,URf) + max(cfL,cfR);

  const scalar SL = min(ULf-cfL, URf-cfR); //Eq. 12 
  const scalar SR = max(ULf+cfL, URf+cfR);

  const scalar pGL = pL+BLsqr/(2.0*mu0);
  const scalar pGR = pR+BRsqr/(2.0*mu0);  

  const scalar SM = (rhoR*URf*(SR-URf) - rhoL*ULf*(SL-ULf) + pGL - pGR)/
    (rhoR*(SR-URf)-rhoL*(SL-ULf));//Eq.38



  /*
    The following implementation is similar to variation HLLC-L (only two-states)
    With the conistent formulation of Mignone (see Jets from Young Stars pp.84)
  */

  //Magnetic Flux
  Bf = (BLf + BRf)/2.0; 
  /*
    It is necessary for the pB equation      
  */


  //Magnetic Flux
  const scalar Bfstar = (SR*BRf-SL*BLf)/(SR-SL); //Eq. 33 (Li), for multi-D


  //------Compute pressure in star region from the right and left sides-------------------
  const scalar pGstarR = pGR + rhoR*(SR-URf)*(SM-URf) - sqr(Bf)+sqr(Bfstar);//Eq.17 Li2005
  const scalar pGstarL = pGL + rhoL*(SL-ULf)*(SM-ULf) - sqr(Bf)+sqr(Bfstar);//Eq.17 Li2005

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




  const scalar pB = (pBR+pBL)/2.0;

  const vector Bhll = (SR*BR-SL*BL + (ULf*BL - Bf*UL + Sf_unit*pB) 
		       - (URf*BR - Bf*UR + Sf_unit*pB))/(SR-SL);

  const scalar rhohll = (SR*rhoR-SL*rhoL + (ULf*rhoL)        
			 - (URf*rhoR))/(SR-SL);

  const vector rhoUhll = (SR*rhoUR-SL*rhoUL + (ULf*rhoUL + pGL*Sf_unit - Bf*BL/mu0) 
			  - (URf*rhoUR + pGR*Sf_unit - Bf*BR/mu0))/(SR-SL);

  /*
    The following way of computing Uhll works better in the MHD vortex case (CFL08)
    However its is worst in the:
    1) 1a shocktube - wingles for CFL08 and 05 
    2) diverge faster in cylinder for CFL05
  */

  const vector Uhll = rhoUhll/rhohll;


  /*
    The following way of computing Uhll works better in the 1D shocktube case
    However its is worst in the MHD vortex at CFL 0.8
  */
  /*
    const vector FBhll = 1.0/(SR-SL)*(
    SR*(ULf*BL - BLf*UL)
    -SL*(URf*BR - BRf*UR)
    +SR*SL*(BR-BL)
    );

    const vector Uhll = (Bhll*SM-FBhll+pB*Sf_unit)/stabilise(Bfstar,SMALL); 
    //Eq. 43 do Mignone (J. f. y. stars)
    //--------------------------------------------------------
    */


  const scalar rhoLstar = rhoL*(SL-ULf)/(SL-SM);
  const scalar rhoRstar = rhoR*(SR-URf)/(SR-SM);

  const vector rhoULstar = 1.0/(SL-SM)*(
					rhoUL*(SL-ULf)
					-(Bfstar*Bhll-Bf*BL)
					-Sf_unit*(pGL-pGstar)
					);

  const vector rhoURstar = 1.0/(SR-SM)*(
					rhoUR*(SR-URf)
					-(Bfstar*Bhll-Bf*BR)
					-Sf_unit*(pGR-pGstar)
					);


  const scalar rhoELstar = (
			    rhoEL*(SL-ULf) - pGL*ULf + pGstar*SM 
			    + Bf*(UL&BL)/mu0 
			    - Bfstar*(Uhll&Bhll)/mu0
			    )/(SL-SM);
                       
  const scalar rhoERstar = (
			    rhoER*(SR-URf) - pGR*URf + pGstar*SM 
			    + Bf*(UR&BR)/mu0 
			    - Bfstar*(Uhll&Bhll)/mu0
			    )/(SR-SM);



  if ( pos(SL) )
    {
      // compute FL
      rhoflux  = (ULf*rhoL)*magSf;
      rhoUflux = (ULf*rhoUL + pGL*Sf_unit - BLf*BL/mu0)*magSf;
      rhoBflux = (ULf*BL - BLf*UL + Sf_unit*pB)*magSf;
      rhoEflux = (ULf*(rhoEL+pGL) - BLf*((UL & BL))/mu0)*magSf;
    }
  else if ( pos(SM) )
    {
      // compute FL*
      rhoflux  = (SM*rhoLstar)*magSf;
      rhoUflux = (rhoULstar*SM + pGstar*Sf_unit - Bfstar*Bhll/mu0)*magSf;
      rhoBflux = (Bhll*SM - Bfstar*Uhll + Sf_unit*pB)*magSf;            
      rhoEflux = ((rhoELstar + pGstar)*SM - Bfstar*(Uhll&Bhll)/mu0)*magSf;
    }
  else if ( pos(SR) )
    {
      // compute FR*
      rhoflux  = (SM*rhoRstar)*magSf;
      rhoUflux = (rhoURstar*SM + pGstar*Sf_unit - Bfstar*Bhll/mu0)*magSf;
      rhoBflux = (Bhll*SM - Bfstar*Uhll + Sf_unit*pB)*magSf;            
      rhoEflux = ((rhoERstar + pGstar)*SM - Bfstar*(Uhll&Bhll)/mu0)*magSf;
    }
  else if ( neg(SR) )
    {
      // compute FR
      rhoflux  = (URf*rhoR)*magSf;
      rhoUflux = (URf*rhoUR + pGR*Sf_unit - BRf*BR/mu0)*magSf;
      rhoBflux = (URf*BR - BRf*UR + Sf_unit*pB)*magSf;
      rhoEflux = (URf*(rhoER+pGR) - BRf*((UR & BR))/mu0)*magSf;
    }
  else
    {
      Info << "Error in HLLD Riemann solver" << endl;
    }




 

  cfmax = max(mag(ULf)+cfL, mag(URf)+cfR)*magSf;

}
