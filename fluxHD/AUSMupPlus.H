/*
  =============================================================================
  AUSM+-up FLUX START

  This subroutine computes the AUSMup+ flux for the Euler equations
  in the direction, S=[Sx,Sy,Sz].
  Conservative form of the Euler equations:

  dU/dt + dF/dx + dG/dy + dH/dz = 0


  Fn = F*Sx + G*Sy + H*Sz = | rho*Uf          |
                            | rho*Uf*u + p*Sx |
                            | rho*Uf*v + p*Sy |
                            | rho*Uf*w + p*Sz |
                            | rho*Uf*H        |    (Uf = u*Sx + v*Sy + w*Sz)


  The mass flux function is given by:

                   |rhoL if Mf>0
  F(i+1/2) = af*Mf<
                   |rhoR otherwise

  where af = min(aL,aR)
  and
  Mf = M4pos + M4neg + Mp.
  M4pos/neg are mach numeber interpolation functions and
  Mp is a pressure diffusion term introduced to enhace calculations at low Mach.

  The pressure flux function is given by:

  pf = P5pos*pL + P5neg*pR + Pu;

  Where P5pos/neg are again Mach number interpolation functions and
  Pu is a velocity diffusion terms

  References
  (1) Meng-Sing Liou, "A sequel to AUSM, Part II: AUSM+-up for all speeds"
  JCP, 2006.

  (2) Xisto, C. M., Páscoa, J. C. , Oliveira, P. J., Nicolini, D. A. (2012),
  An hybrid pressure­ density based algorithm for the Euler equations at
  all Mach number regimes, International Journal for Numerical Methods
  in Fluids Vol. 70:8, pages 961-976, doi:10.1002/fld.2722.


  Authors: Carlos Xisto (xisto@ubi.pt)
  Date: ..-..-2014
  =============================================================================
*/

void AUSMupPlusFlux (
		     const vector &Sf,
		     const scalar &magSf,
		     const scalar &rhoL,  const scalar &rhoR,
		     const vector &rhoUL, const vector &rhoUR,
		     const vector &UL,    const vector &UR,
		     const scalar &pL,    const scalar &pR,
		     const scalar &rhoEL, const scalar &rhoER,
		     const scalar &gammaL, const scalar &gammaR,
		     const scalar &Ms,
		     scalar &rhoflux,   vector &rhoUflux,
		     scalar &rhoEflux, scalar &cfmax
		     )
{


  const scalar Ku = 0.75;
  const scalar Kp = 0.25;
  const scalar beta = 1.0/8.0;


  //================== Direction Spliting Uf = U & Sf    =======================
  const vector Sf_unit = Sf /magSf;

  const scalar    ULf = UL & Sf_unit;
  const scalar    URf = UR & Sf_unit;

  //============================================================================
  //    CALCULATE NUMERICAL SOUND SPEED (Eq(29, liou2006))
  //=============================================================================

  //------------ TOTAL ENTALPHY -------------------------------------------------

  const scalar HtL = rhoEL + pL;
  const scalar HtR = rhoER + pR;

  const scalar Htf = min(HtL/rhoL, HtR/rhoR);

  const scalar gammaF = 0.5*(gammaL+gammaR);


  //------------NUMERICAL SOUND SPEED and  MACH NUMBER Eq(16, liou2006)----------
  const scalar acritF = Foam::sqrt(mag((2.0*Htf*(gammaF-1.0))/(gammaF+1.0)));

  const scalar    aL = sqr(acritF)/(max(acritF, ULf));   //Eq(30, liou2006)
  const scalar    aR = sqr(acritF)/(max(acritF,-URf));

  const scalar   af = min(aL, aR);

  const scalar  ML = ULf/af;
  const scalar  MR = URf/af;

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
    * ((pR - pL)/(rhoM*sqr(af))); //(Eq73)

  const scalar   Pu = -Ku*P5pos*P5neg*(rhoL+rhoR)
    *(fpre*(af))*((URf-ULf));   //termo difusivo (Eq75)


  //=============================================================================
  //                                  FLUX COMPUTATION
  //=============================================================================

  const scalar Uf = af*(M4pos+M4neg+Mp);
  const scalar pf = pL*P5pos + pR*P5neg + Pu;

  rhoflux =  0.5*(Uf*(rhoL+rhoR)-mag(Uf)*(rhoR-rhoL))*magSf;

  rhoUflux = 0.5*(Uf*(rhoUL+rhoUR)-mag(Uf)*(rhoUR-rhoUL))*magSf + pf*Sf;

  rhoEflux = 0.5*(Uf*(HtL+HtR)-mag(Uf)*(HtR-HtL))*magSf;

  //Co variables
  cfmax = (mag(Uf)+af)*magSf;

}
//=============================================================================
//                                   AUSM+ FLUX END
//=============================================================================
