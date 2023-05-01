/*
  =============================================================================
  HLLC FLUX START

  This subroutine computes the HLLC flux for the Euler equations
  in the direction, S=[Sx,Sy,Sz].
  Conservative form of the Euler equations:

  dU/dt + dF/dx + dG/dy + dH/dz = 0

  Fn = F*Sx + G*Sy + H*Sz = | rho*Uf          |
                            | rho*Uf*u + p*Sx |
                            | rho*Uf*v + p*Sy |
                            | rho*Uf*w + p*Sz |
                            | rho*Uf*H        |    (Uf = u*Sx + v*Sy + w*Sz)


  The flux function is given by:


           |Fn(UL)  if SL>0,
           |F(UL*)  if SL<0<SM,
  FHLLC = <
           |F(UR*) if SM<0<SR,
           |Fn(UR)  if Sr<0,


  SL=min(Uf-af, UfTile-aTilde);
  SR=max(Uf+af, UfTile+aTilde);


  References
  (1) Toro, E. F.; Spruce, M.; Speares, W. (1994), "Restoration of the contact
  surface in the HLL-Riemann solver", Shock Waves 4: 25–34,doi:10.1007/BF01414629

  (2)P. Batten, M.A. Leschziner, U.C. Goldberg, Average-State Jacobians and
  Implicit Methods for Compressible Viscous and Turbulent Flows, Journal of
  Computational Physics, Volume 137, Issue 1, October 1997, Pages 38-78,
  ISSN 0021-9991, http://dx.doi.org/10.1006/jcph.1997.5793.


  Authors: Carlos Xisto (xisto@ubi.pt)
  Date:18-02-2015
  =============================================================================
*/


void HLLCFlux (
	       const vector &Sf,
	       const scalar &magSf,
	       const scalar &rhoL,  const scalar &rhoR,
	       const vector &rhoUL, const vector &rhoUR,
	       const vector &UL,    const vector &UR,
	       const scalar &pL,    const scalar &pR,
	       const scalar &rhoEL, const scalar &rhoER,
	       const scalar &gammaL, const scalar &gammaR,
	       scalar &rhoflux,   vector &rhoUflux,
	       scalar &rhoEflux, scalar &cfmax
	       )

{


  //--------------------- Direction Spliting Uf = U & Sf | Bf = B & Sf ---------
  const vector Sf_unit = Sf /magSf;

  const scalar ULf = UL & Sf_unit;
  const scalar URf = UR & Sf_unit;

  //-------FAST MAGNETOSONIC SPEED----------------------------------------------
  const scalar    aL = Foam::sqrt(gammaL*pL/rhoL);
  const scalar    aR = Foam::sqrt(gammaR*pR/rhoR);

  //-------Signal speeds--------------------------------------------------------

  // Compute left and right total enthalpies:
  const scalar HL  = (rhoEL + pL)/rhoL;
  const scalar HR =  (rhoER + pR)/rhoR;

  // compute Roe weights for Eq. 13
  const scalar rhoLSqrt = Foam::sqrt(rhoL);
  const scalar rhoRSqrt = Foam::sqrt(rhoR);

  const scalar wL = rhoLSqrt/(rhoLSqrt + rhoRSqrt);
  const scalar wR = 1.0 - wL;

  // Roe averaged velocity
  const vector UTilde = UL*wL + UR*wR;
  const scalar UfTilde = UTilde & Sf_unit;


  // Roe averaged total enthalpy
  const scalar HTilde = wL*HL + wR*HR;

  // Roe averaged gamma
  //const scalar gammaTilde = wL*gammaL + wR*gammaR;
  const scalar gammaTilde = gammaL;


  // Speed of sound with Roe reconstruction values
  const scalar aTilde =
    Foam::sqrt((gammaTilde-1.0)*(HTilde-0.5*magSqr(UTilde)));

  const scalar SL = min(ULf-aL, UfTilde-aTilde); //Eq. 13
  const scalar SR = max(URf+aR, UfTilde+aTilde);

  //const scalar SL = min(ULf-aL, URf-aR); //Eq. 12
  //const scalar SR = max(ULf+aL, URf+aR);

  //stabilise((rhoR*(SR-URf)-rhoL*(SL-ULf)),SMALL); deve ajudar a resolver sigfe
  const scalar SM = (rhoR*URf*(SR-URf)-rhoL*ULf*(SL - ULf) + pL - pR)/
    (rhoR*(SR-URf)-rhoL*(SL-ULf));//Eq.14



  //------Compute pressure in star region from the right side-------------------
  const scalar pStarR = rhoR * (URf - SR) * (URf - SM) + pR;//Eq.16

  //------Should be equal to the left and intermediate sides--------------------
  const scalar pStarL = rhoL * (ULf - SL) * (ULf - SM) + pL; //Eq.16

  // give a warning if this is not the case
  if ( mag(pStarR-pStarL) > 1e-6 )
    {
      Info << "WARNING, mag(pStarRight-pStarLeft) > SMALL " << endl;
    }

  //use pStarRight for pStar, as in theory, pTL*=PTL**=pTR**=pTR*=pT*
  const scalar pStar = pStarR;


  /*
    The following implementation is similar to variation 2 of HLLC
    http://www.prague-sum.com/download/2012/Toro_2-HLLC-RiemannSolver.pdf

    Fk*=(S* (Sk Uk - Fk) + Sk p* D*)/(Sk - S*),

    where k=L,R and p* is the average between pL* and pR*
    However in principle pL*=pR*
  */

  if ( pos(SL) )
    {
      // compute FL
      rhoflux  = (ULf*rhoL)*magSf;
      rhoUflux = (ULf*rhoUL+pL*Sf_unit)*magSf;
      rhoEflux = (ULf*(rhoEL+pL))*magSf;
    }
  else if ( pos(SM) )
    {
      const scalar omegaL=1.0/(SL-SM);
      // compute FL*
      rhoflux  = (
		  omegaL*(SM*rhoL*(SL-ULf))
		   )*magSf;

      rhoUflux = (
		  omegaL*(SM*((SL-ULf)*rhoUL
			      +(pStar-pL)*Sf_unit))+pStar*Sf_unit
		   )*magSf;

      rhoEflux = (
		  omegaL *(SM*((SL-ULf)*rhoEL-pL*ULf
			       +pStar*SM)) + pStar*SM
		  )*magSf;

    }
  else if ( pos(SR) )
    {
      const scalar omegaR=1.0/(SR-SM);
      // compute FR*
      rhoflux  = (
		  omegaR*(SM*rhoR*(SR-URf))
		   )*magSf;
      rhoUflux = (
		  omegaR*(SM*((SR-URf)*rhoUR
			      +(pStar-pR)*Sf_unit))+pStar*Sf_unit
		   )*magSf;
      rhoEflux = (
		  omegaR *(SM*((SR-URf)*rhoER-pR*URf
			       +pStar*SM)) + pStar*SM
		  )*magSf;
    }
  else if ( neg(SR) )
    {
      // compute FR

      rhoflux  = (URf*rhoR)*magSf;
      rhoUflux = (URf*rhoUR+pR*Sf_unit)*magSf;
      rhoEflux = (URf*(rhoER+pR))*magSf;

    }
  else
    {
      Info << "Error in HLLC Riemann solver" << endl;
    }


  cfmax = max(mag(ULf)+aL, mag(URf)+aR)*magSf;
}
