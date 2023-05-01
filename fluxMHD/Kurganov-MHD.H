/*
  =============================================================================
  Kurganov FLUX START

  This subroutine computes the KT flux for the MHD equations
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



  Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
  flux description. Therefore do not stress to much if you find a pB*Sf below.

  References
  (1) A. Kurganov, E. Tadmor, New High-Resolution Central Schemes for Nonlinear
  Conservation Laws and Convection–Diffusion Equations, Journal of Computational
  Physics, Volume 160, Issue 1, 1 May 2000, Pages 241-282, ISSN 0021-9991,
  http://dx.doi.org/10.1006/jcph.2000.6459.


  Modification for MHD:
  Replacment of acoutic wave:

  |Uf|+af

  by the magnetosonic fast wave:

  |Uf|+cf

  Authors: Carlos Xisto (xisto@ubi.pt)
  Date: 13-02-2015
  =============================================================================
*/
void KurganovMHDFlux (
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

  //-------------------Direction Spliting Uf = U & Sf--------------------------------
  const vector Sf_unit = Sf/magSf;

  const scalar ULf = UL & Sf_unit;
  const scalar URf = UR & Sf_unit;

  const scalar BLf = BL  & Sf_unit;
  const scalar BRf = BR  & Sf_unit;



  //-------FAST MAGNETOSONIC SPEED-------------------------------------------------------
  const scalar    aL = Foam::sqrt(gammaL*pL/rhoL);
  const scalar    aR = Foam::sqrt(gammaR*pR/rhoR);

  const scalar BLsqr = magSqr(BL);
  const scalar BRsqr = magSqr(BR);

  const scalar bfL = BLsqr/(rhoL*mu0); //alfven speed
  const scalar bfR = BRsqr/(rhoR*mu0); //alfven speed

  const scalar caL = BLf/Foam::sqrt(rhoL*mu0); //normal alfven speed
  const scalar caR = BRf/Foam::sqrt(rhoR*mu0);


  const scalar cfL = Foam::sqrt((1.0/2.0) * (sqr(aL) + bfL + Foam::sqrt(sqr(sqr(aL)+ bfL)
									- 4.0*sqr(aL)*sqr(caL) )));//left value for magnetosonic fast speed
  const scalar cfR = Foam::sqrt((1.0/2.0) * (sqr(aR) + bfR + Foam::sqrt(sqr(sqr(aR)+ bfR)
									- 4.0*sqr(aR)*sqr(caR) )));//rigth value for magnetosonic fast speed


  const scalar ap = max(max(ULf + cfL, URf + cfR), 0.0);
  const scalar am = min(min(ULf - cfL, URf - cfR), 0.0);
  const scalar aPos=(ap/(ap - am));
  //const scalar amaxSf=max(mag(am), mag(ap)); it will be necessary for tadmor
  const scalar aSf=(am*aPos);


  const scalar aNeg=(1.0 - aPos);

  const scalar ULfn = ULf*aPos;
  const scalar URfn = URf*aNeg;

  const scalar cL=(ULfn - aSf);
  const scalar cR=(URfn + aSf);


  const scalar amaxSf = max(mag(cL), mag(cR));

  //Magnetic Flux
  Bf = (BLf + BRf)/2.0;
  /*
    It is necessary for the pB equation
  */


  //CONTINUITY FLUX
  rhoflux = (cL*rhoL + cR*rhoR)*magSf;


  //MOMENTUM FLUX
  const scalar pGL = pL+BLsqr/(2.0*mu0);
  const scalar pGR = pR+BRsqr/(2.0*mu0);

  const vector BfB = (BLf*BL +  BRf*BR)/(2.0*mu0);

  rhoUflux = ( (cL*rhoUL + cR*rhoUR)
	       + (aPos*pGL + aNeg*pGR)*Sf_unit - BfB)*magSf;


  //INDUCTION FLUX
  const scalar  pBf = 0.5*(pBL + pBR);
  const vector  UfB = (cL*BL + cR*BR);
  const vector  BfU = (aPos*( BLf*UL)+aNeg*( BRf*UR));
  rhoBflux = (UfB-BfU)*magSf + pBf*Sf;


  //ENERGY FLUX

  const scalar BUL = UL & BL;
  const scalar BUR = UR & BR;
  const scalar BfBU = aPos*( BLf*BUL)/mu0 + aNeg*( BRf*BUR)/mu0;

  rhoEflux = (cL*(rhoEL + pGL) + cR*(rhoER + pGR)
	      + aSf*pGL - aSf*pGR - BfBU)*magSf;



  cfmax =amaxSf*magSf;

  //============================================================================
  //                   kurganovMHD END
  //============================================================================

}
