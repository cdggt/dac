/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DirACFoamWallReduce.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "fvcSmooth.H"
#include "fvcSurfaceIntegrate.H"

#include "fvcDdt.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamWallReduce::thermophysicalPredictor()
{

    //Make references to all volume field variables
    volScalarField& Xl(Xl_);
    volScalarField& Xv(Xv_);
    volScalarField& T(T_);
    volScalarField& j(j_);

    //Hard coded conversion factors/unit stuff
    dimensionedScalar AP("AP", dimless, 101325);
    dimensionedScalar CU("CU", dimTemperature, 1);
    dimensionedScalar jConst("jConst", dimless/dimTime, 0);

    //Calculations necessary for the coefficients in the transport equations
    volScalarField Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    volScalarField Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    volScalarField dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    volScalarField ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    Xv = (1/AP)*pow(10,A-((B*CU)/((T/(1-(dP/(rhoL*Lw))))+(C*CU))));

    volScalarField dXvdT = ((Xv*rhoL*Lw*std::log(10)*pow(A-log10(AP*Xv),2))/(CU*B*((rhoL*Lw)-(dP))));

    volScalarField dXvdXl = -dXvdT*(T/(2*Xi*(1-Xs-Xa)*(rhoL*Lw-dP)))*ddPdXi;

    volScalarField Dfact = (1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))));

    volScalarField KappaW = Xs*kappaS + Xa*kappaA + Xl*kappaL;

    volScalarField cW = Xs*CpS*rhoS+Xa*CpA*rhoA+Xl*CpL*rhoL;

    volScalarField psi = rhoV*(1 - Xs - Xa - Xl)*Lw;

    volScalarField gam = (rhoV/rhoL)*(1 - Xs - Xa - Xl);

    j=dXvdT*fvc::ddt(T)-fvc::laplacian((dXvdT/(tau*tau))*Dfact*D0,T)+dXvdXl*fvc::ddt(Xl)-fvc::laplacian(dXvdXl*Dfact*D0,Xl)/(tau*tau);



    //Solve for liquid transport
    fvScalarMatrix XlEqn =
    (
        (1+(gam*dXvdXl))*fvm::ddt(Xl)
      + gam*dXvdT*fvc::ddt(T)
      ==
        fvm::laplacian((((((1-Xs-Xa)*pow(r0,2))/(8*eta))*Co*ddPdXi*(-1/(2*(1-Xs-Xa)*Xi)))+(gam*dXvdXl*Dfact*D0))/(tau*tau),Xl)
      + fvc::laplacian(gam*dXvdT*Dfact*D0/(tau*tau),T)
    );

    XlEqn.relax();

    fvConstraints().constrain(XlEqn);

    XlEqn.solve("Xl");



    Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    Xv = (1/AP)*pow(10,A-((B*CU)/((T/(1-(dP/(rhoL*Lw))))+(C*CU))));

    dXvdT = ((Xv*rhoL*Lw*std::log(10)*pow(A-log10(AP*Xv),2))/(CU*B*((rhoL*Lw)-(dP))));

    dXvdXl = -dXvdT*(T/(2*Xi*(1-Xs-Xa)*(rhoL*Lw-dP)))*ddPdXi;

    Dfact = (1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))));

    KappaW = Xs*kappaS + Xa*kappaA + Xl*kappaL;

    cW = Xs*CpS*rhoS+Xa*CpA*rhoA+Xl*CpL*rhoL;

    psi = rhoV*(1 - Xs - Xa - Xl)*Lw;

    gam = (rhoV/rhoL)*(1 - Xs - Xa - Xl);

    j=dXvdT*fvc::ddt(T)-fvc::laplacian((dXvdT/(tau*tau))*Dfact*D0,T)+dXvdXl*fvc::ddt(Xl)-fvc::laplacian(dXvdXl*Dfact*D0,Xl)/(tau*tau);



    //Solve for heat transport
    fvScalarMatrix EEqn =
    (
        (cW+psi*dXvdT)*fvm::ddt(T)
      + ((CpL*rhoL*T)+(psi*dXvdXl))*fvc::ddt(Xl)
      ==
        fvm::laplacian(KappaW+((psi*dXvdT/(tau*tau))*Dfact*D0),T)
      + fvc::laplacian(((psi*dXvdXl/(tau*tau))*Dfact*D0),Xl)
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();



    Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    Xv = (1/AP)*pow(10,A-((B*CU)/((T/(1-(dP/(rhoL*Lw))))+(C*CU))));

    dXvdT = ((Xv*rhoL*Lw*std::log(10)*pow(A-log10(AP*Xv),2))/(CU*B*((rhoL*Lw)-(dP))));

    dXvdXl = -dXvdT*(T/(2*Xi*(1-Xs-Xa)*(rhoL*Lw-dP)))*ddPdXi;

    Dfact = (1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))));

    KappaW = Xs*kappaS + Xa*kappaA + Xl*kappaL;

    cW = Xs*CpS*rhoS+Xa*CpA*rhoA+Xl*CpL*rhoL;

    psi = rhoV*(1 - Xs - Xa - Xl)*Lw;

    gam = (rhoV/rhoL)*(1 - Xs - Xa - Xl);

    j=dXvdT*fvc::ddt(T)-fvc::laplacian((dXvdT/(tau*tau))*Dfact*D0,T)+dXvdXl*fvc::ddt(Xl)-fvc::laplacian(dXvdXl*Dfact*D0,Xl)/(tau*tau);



    volScalarField diffXl = Xl-Xl.oldTime();

    dXlAvg = gAverage(diffXl);
    dXlNum = max(gMax(diffXl),-gMin(diffXl));

    Info<< "Xl Difference mean: " << dXlAvg << " max: " << dXlNum << endl;



    volScalarField diffT = T-T.oldTime();

    dTAvg = gAverage(diffT);
    dTNum = max(gMax(diffT),-gMin(diffT));

    Info<< "T Difference mean: " << dTAvg << " max: " << dTNum << endl;
}


// ************************************************************************* //
