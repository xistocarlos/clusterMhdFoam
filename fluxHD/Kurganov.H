/*
=============================================================================
Kurganov FLUX START

 This subroutine computes the KT flux for the Euler equations
 in the direction, S=[Sx,Sy,Sz].
 Conservative form of the Euler equations:

     dU/dt + dF/dx + dG/dy + dH/dz = 0


     Fn = F*Sx + G*Sy + H*Sz = | rho*Uf          |
                               | rho*Uf*u + p*Sx |
                               | rho*Uf*v + p*Sy |
                               | rho*Uf*w + p*Sz |
                               | rho*Uf*H        |    (Uf = u*Sx + v*Sy + w*Sz)

The flux function is given by:





References
(1) A. Kurganov, E. Tadmor, New High-Resolution Central Schemes for Nonlinear
    Conservation Laws and Convection–Diffusion Equations, Journal of Computational
    Physics, Volume 160, Issue 1, 1 May 2000, Pages 241-282, ISSN 0021-9991,
    http://dx.doi.org/10.1006/jcph.2000.6459.

Authors: Carlos Xisto (xisto@ubi.pt)
Date:1-02-2015
=============================================================================
*/
void KurganovFlux (
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



//-------------------Direction Spliting Uf = U & Sf--------------------------------
const vector Sf_unit = Sf/magSf;

const scalar ULf = UL & Sf_unit;
const scalar URf = UR & Sf_unit;


//speed of sound
const scalar    aL = Foam::sqrt(gammaL*pL/rhoL);
const scalar    aR = Foam::sqrt(gammaR*pR/rhoR);


const scalar ap = max(max(ULf + aL, URf + aR), 0.0);
const scalar am = min(min(ULf - aL, URf - aR), 0.0);
const scalar aPos=(ap/(ap - am));
//const scalar amaxSf=max(mag(am), mag(ap)); it will be necessary for tadmor
const scalar aSf=(am*aPos);


const scalar aNeg=(1.0 - aPos);

const scalar ULfn = ULf*aPos;
const scalar  URfn = URf*aNeg;

const scalar cL=(ULfn - aSf);
const scalar cR=(URfn + aSf);


const scalar amaxSf = max(mag(cL), mag(cR));


//CONTINUITY FLUX
rhoflux = (cL*rhoL + cR*rhoR)*magSf;


//MOMENTUM FLUX
 rhoUflux = ( (cL*rhoUL + cR*rhoUR)
          + (aPos*pL + aNeg*pR)*Sf_unit )*magSf;


//ENERGY FLUX
  rhoEflux = (cL*(rhoEL + pL) + cR*(rhoER + pR)
          + aSf*pL - aSf*pR)*magSf;



  cfmax =amaxSf*magSf;

}
