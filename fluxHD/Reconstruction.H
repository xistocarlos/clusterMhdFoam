/*
=============================================================================
                LEFT & RIGTH FACE VALUES

Using J. Greenshields approach

Christopher J. Greenshields, Henry G. Weller, Luca Gasparini3 and Jason M. Reese1.
"Implementation of semi-discrete, non-staggered central schemes in a colocated,
polyhedral, finite volume framework,for high-speed viscous flows", IJNMF,2009.
==================================================================================
*/

 surfaceVectorField    rhoUL = fvc::interpolate(rhoU, pos, "reconstruct(U)");
 surfaceVectorField    rhoUR = fvc::interpolate(rhoU, neg, "reconstruct(U)");

 surfaceScalarField    rhoL = fvc::interpolate(rho, pos, "reconstruct(rho)");
 surfaceScalarField    rhoR = fvc::interpolate(rho, neg, "reconstruct(rho)");

 surfaceScalarField    pL = fvc::interpolate(p, pos, "reconstruct(p)");
 surfaceScalarField    pR = fvc::interpolate(p, neg, "reconstruct(p)");

 surfaceScalarField    rhoEL = fvc::interpolate(rhoE, pos, "reconstruct(T)");
 surfaceScalarField    rhoER = fvc::interpolate(rhoE, neg, "reconstruct(T)");

 volScalarField gamma = thermo.Cp()/thermo.Cv();

 surfaceScalarField    gammaL = fvc::interpolate(gamma, pos, "reconstruct(T)");
 surfaceScalarField    gammaR = fvc::interpolate(gamma, neg, "reconstruct(T)");


/*
=============================================================================
//Linear Recontruction with slope limiter
//To BE DONE
==============================================================================
*/
