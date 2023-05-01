word fluxScheme("");
if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{
  if
  (
    (fluxScheme == "TVDLF-MHD") ||
    (fluxScheme == "AUSMpw-MHD") ||
    (fluxScheme == "AUSMup-MHD") ||
    (fluxScheme == "Kurganov-MHD") ||
    (fluxScheme == "HLLD") ||
    (fluxScheme == "HLLC-L") ||
    (fluxScheme == "HLLC") ||
    (fluxScheme == "Roe")  ||
    (fluxScheme == "AUSMupPlus") ||
    (fluxScheme == "TVDLF") ||
    (fluxScheme == "Kurganov") ||
    (fluxScheme == "Roe-MHD")


  )

  {
    Info<< "fluxScheme is: " << fluxScheme << nl
    << " HLLC, Roe, AUSMupPlus, TVDLF, and Kurganov not working with MHD"
    << nl
    << endl;
  }
  else
  {
    FatalErrorIn
    (
      "clusterFoam::readFluxScheme"
    )   << "fluxScheme: " << fluxScheme
    << " is not a valid choice. " << nl
    << "Options are: AUSMpw-MHD; TVDLF-MHD;" << nl
    << "             HLLC(GD); Roe-MHD(NI); AUSMup-MHD"<< nl
    << "             Kurganov-MHD; HLLC-L; HLLD; Roe(GD)"<< nl
    << "             Kurganov(GD); AUSMupPlus(GD); TVDLF(GD)"<< nl
    << abort(FatalError);
  }
}
