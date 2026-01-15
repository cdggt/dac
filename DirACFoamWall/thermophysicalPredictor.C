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

#include "DirACFoamWall.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "fvcSmooth.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcLaplacian.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamWall::thermophysicalPredictor()
{

    //Make references to all volume field variables
    volScalarField& Xl(Xl_);
    volScalarField& Xv(Xv_);
    volScalarField& T(T_);
    volScalarField& j(j_);

    //Hard coded conversion factors/unit stuff
    dimensionedScalar AP("AP", dimless, 101325);
    dimensionedScalar CU("CU", dimTemperature, 1);

    //Calculations necessary for the coefficients in the transport equations
    volScalarField Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    volScalarField Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    volScalarField dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    volScalarField ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    volScalarField Ts = CU*((B/(A-log10(max(1e-9,Xv*AP)))) - C);

    j = fact*(2/(r0*Xi))*sqrt((Rw*max(T,CU))/(2*Foam::constant::mathematical::pi))
        *((dP/(rhoL*Rw*max(T,CU)))+((Lw/Rw)*((1/Ts)-(1/max(T,CU)))));



    //Solve for liquid transport
    fvScalarMatrix XlEqn =
    (
        fvm::ddt(Xl)
      + fvm::div(phi,Xl)
      ==
        fvm::laplacian((((1-Xs-Xa)*pow(r0,2))/(8*eta))*Co*ddPdXi*(-1/(2*(1-Xs-Xa)*Xi))/(tau*tau),Xl)
      - (1 - Xs - Xa - Xl)*(rhoV/rhoL)*j
    );

    XlEqn.relax();

    fvConstraints().constrain(XlEqn);

    XlEqn.solve("Xl");

    fvConstraints().constrain(Xl);



    Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    Ts = CU*((B/(A-log10(max(1e-9,Xv*AP)))) - C);

    j = fact*(2/(r0*Xi))*sqrt((Rw*max(T,CU))/(2*Foam::constant::mathematical::pi))
      *((dP/(rhoL*Rw*max(T,CU)))+((Lw/Rw)*((1/Ts)-(1/max(T,CU)))));



    //Solve for vapor transport
    fvScalarMatrix XvEqn =
    (
        fvm::ddt(Xv)
      + fvm::div(phi,Xv)
      ==
        fvm::laplacian((1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))))*D0/(tau*tau),Xv)
      + j
    );

    XvEqn.relax();

    fvConstraints().constrain(XvEqn);

    XvEqn.solve("Xv");

    fvConstraints().constrain(Xv);

    

    Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    Ts = CU*((B/(A-log10(max(1e-9,Xv*AP)))) - C);

    j = fact*(2/(r0*Xi))*sqrt((Rw*max(T,CU))/(2*Foam::constant::mathematical::pi))
      *((dP/(rhoL*Rw*max(T,CU)))+((Lw/Rw)*((1/Ts)-(1/max(T,CU)))));



    //Solve for heat transport
    fvScalarMatrix EEqn =
    (
       (Xs*CpS*rhoS+Xa*CpA*rhoA+Xl*CpL*rhoL)*fvm::ddt(T)
     + (CpL*rhoL*T)*fvc::ddt(Xl)
     + (Xs*CpS*rhoS+Xa*CpA*rhoA+Xl*CpL*rhoL)*fvm::div(phi,T)
     ==
       fvm::laplacian(Xs*kappaS + Xa*kappaA + Xl*kappaL,T)
     - rhoV*(1 - Xs - Xa - Xl)*Lw*j
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(T);



    //Recalculate j for use in time step control (the benefit of more precise time
    //stepping outweighs the small extra computational cost)
    Xi = sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa));

    Co = 1 - 4*pow(Xi,2) + 3*pow(Xi,4) - 4*pow(Xi,4)*log(Xi);

    dP = -((sig)/(r0*Xi)) - ((AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,3)));

    ddPdXi = ((sig)/(r0*pow(Xi,2))) - ((AH)/(2*Foam::constant::mathematical::pi*pow(r0,3)*pow(1-Xi,4)));

    Ts = CU*((B/(A-log10(max(1e-9,Xv*AP)))) - C);

    j = fact*(2/(r0*Xi))*sqrt((Rw*max(T,CU))/(2*Foam::constant::mathematical::pi))
      *((dP/(rhoL*Rw*max(T,CU)))+((Lw/Rw)*((1/Ts)-(1/max(T,CU)))));



    volScalarField diffXv = Xv-Xv.oldTime();

    dXvAvg = gAverage(diffXv);
    dXvNum = max(gMax(diffXv),-gMin(diffXv));

    Info<< "Xv Difference mean: " << dXvAvg << " max: " << dXvNum << endl;



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
