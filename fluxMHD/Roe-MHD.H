/*
  =============================================================================
  ROE FLUX START

  This subroutine computes the Roe flux for the MHD equations
  in the direction, S=[Sx,Sy,Sz].
  Conservative form of the Euler equations:

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

  F(i+1/2)  = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ],

  where

  An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.

  The dissipation term, |An|dU, is actually computed as

  sum_{k=1,5} |lambda_k| * (LdU)_k * r_k,

  where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
  and r_k is the k-th right-eigenvector.

  Warning: The Dedner hyperbolic/parabolic correction for div(B) is not present in the
  flux description. Therefore do not stress to much if you find a pB*Sf below.

  References:
  (1)


  NOT IMPLEMENTED YET..
  =============================================================================
*/


void RoeMHDFlux (
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
		 scalar &Bflux, scalar &rhoflux,
		 vector &rhoUflux, vector &rhoBflux,
		 scalar &rhoEflux, scalar &cfmax
		 )

{

  //--------------------- Direction Spliting Uf = U & Sf | Bf = B & Sf ---------
  // Direction Spliting Uf = U & Sf
  const vector Sf_unit = Sf/magSf;

  const scalar ULf = UL & Sf_unit;
  const scalar URf = UR & Sf_unit;


  // Compute left and right total enthalpies:
  const scalar HtL  = (rhoEL + pL)/rhoL;
  const scalar HtR =  (rhoER + pR)/rhoR;

  // ROE averaged quantities
  const scalar RT = Foam::sqrt(rhoR/rhoL);

  const scalar   rhoTilde = RT*rhoL;                     //Roe-averaged density
  const vector     UTilde = (UL + RT*UR)/(1.0 + RT);     //Roe-averaged Velocity
  const scalar    UfTilde = (UTilde & Sf_unit);          //Roe-averaged face-normal veloc
  const scalar     HTilde = (HtL + RT*HtR)/(1.0 + RT);   //Roe-averaged total enthalpy
  const scalar gammaTilde = gammaL;                      //Roe-averaged gamma
  const scalar aTilde = Foam::sqrt((gammaTilde - 1.0)    //Roe-averaged speed of sound
				   *(HTilde-0.5*magSqr(UTilde)));

  // Primitive diferences
  const scalar deltaRho = rhoR - rhoL;        //Density difference
  const scalar deltaP   =   pR - pL;          //Pressure difference
  const vector deltaU   =   UR - UL;          //Velocity difference
  const scalar deltaUf  = (deltaU & Sf_unit); //Normal Velocity difference


  //Calculate wave Strengths (alpha)
  const scalar alpha1 = (deltaP - aTilde*rhoTilde*deltaUf)
    /(2.0*sqr(aTilde));                //Left-moving acoustic wave strength
  const scalar alpha2 = deltaRho - deltaP/sqr(aTilde);     //Entropy wave strength
  const scalar alpha3 = (deltaP + aTilde*rhoTilde*deltaUf)
    /(2.0*sqr(aTilde));                //Right-moving acoustic wave strength
  const scalar alpha4 = rhoTilde;    //Shear wave strength (not really, just a factor)


  //Roe Matrix eingenvalues
  scalar l1 = mag(UfTilde - aTilde); //Left-moving acoustic wave
  scalar l2 = mag(UfTilde);          //Entropy wave
  scalar l3 = mag(UfTilde + aTilde); //Right-moving acoustic wave
  scalar l4 = mag(UfTilde);          //Shear waves

  /*
    Harten's Entropy Fix JCP(1983), 49, pp357-393.
  */
  const scalar eps = 0.5*aTilde; //com 0.5*aTilde funciona no blunt (o valor poderá ser mais baixo)

  if(l1 < eps)
    {
      l1 = (sqr(l1) + sqr(eps))/(2.0*eps);
    }

  if(l2 < eps)
    {
      l2 = (sqr(l2) + sqr(eps))/(2.0*eps);
    }

  if(l4 < eps)
    {
      l4 = (sqr(l4) + sqr(eps))/(2.0*eps);
    }

  if(l3 < eps)
    {
      l3 = (sqr(l3) + sqr(eps))/(2.0*eps);
    }

  //Right Eigenvectors

  //Left-moving acoustic wave
  const scalar R11 = 1.0;
  const vector R21 = UTilde - aTilde*Sf_unit;
  const scalar R31 = HTilde - aTilde*UfTilde;

  //Entropy wave
  const scalar R12 = 1.0;
  const vector R22 = UTilde;
  const scalar R32 = 0.5*magSqr(UTilde);


  //Rigth-moving acoustic wave
  const scalar R13 = 1.0;
  const vector R23 = UTilde + aTilde*Sf_unit;
  const scalar R33 = HTilde + aTilde*UfTilde;

  //Shear wave
  const scalar R14 = 0.0;
  const vector R24 = deltaU - deltaUf*Sf_unit;
  const scalar R34 = (UTilde & deltaU) - deltaUf*UfTilde;


  //Dissipation Terms:
  //|An|(UR-UL) = R|Lambda|L*dU = sum_k of [l(k)*alpha(k)* R(:,k)]
  const scalar  dF11 = l1*alpha1*R11;
  const vector dF124 = l1*alpha1*R21;
  const scalar  dF15 = l1*alpha1*R31;

  const scalar  dF21 = l2*alpha2*R12;
  const vector dF224 = l2*alpha2*R22;
  const scalar  dF25 = l2*alpha2*R32;

  const scalar  dF31 = l3*alpha3*R13;
  const vector dF324 = l3*alpha3*R23;
  const scalar  dF35 = l3*alpha3*R33;

  const scalar  dF41 = l4*alpha4*R14;
  const vector dF424 = l4*alpha4*R24;
  const scalar  dF45 = l4*alpha4*R34;

  //Left Flux
  const scalar   phiL = rhoL*ULf;
  const vector phiUpL = rhoUL*ULf + pL*Sf_unit;
  const scalar phiEpL = (rhoEL + pL)*ULf;

  //Rigth Flux
  const scalar   phiR = rhoR*URf;
  const vector phiUpR = rhoUR*URf + pR*Sf_unit;
  const scalar phiEpR = (rhoER + pR)*URf;


  //Compute Final Fluxes
  rhoflux = 0.5*(phiL + phiR - (dF11 + dF21 + dF31 + dF41))*magSf;
  rhoUflux = 0.5*(phiUpL + phiUpR - (dF124 + dF224 + dF324 + dF424))*magSf;
  rhoEflux = 0.5*(phiEpL + phiEpR - (dF15 + dF25 + dF35 + dF45))*magSf;

  //=============================================================================
  //ROE FLUX END
  //=============================================================================


  // CFL variables
  cfmax = (mag(UfTilde) + aTilde)*magSf;
}
