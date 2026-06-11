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

#include "DirACFoamWallAdsorb.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "fvcSmooth.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamWallAdsorb::thermophysicalPredictor()
{

    //Make references to all volume field variables
    volScalarField& CO2(CO2_);
    volScalarField& CO2a(CO2a_);

    volScalarField Zero = CO2*0;
    dimensionedScalar tFix("tFix", dimTime, 1);


    /*
    //Solve for vapor transport
    fvScalarMatrix CO2aEqn =
    (
        fvm::ddt(CO2a)
      + fvm::div(phi,CO2a)
      ==
        (alpha*neg(CO2a-cap))*fvc::laplacian(((1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))))*D0/(3*tau)),CO2)
    );

    CO2aEqn.relax();

    fvConstraints().constrain(CO2aEqn);

    CO2aEqn.solve("CO2a");

    fvConstraints().constrain(CO2a);
    */


    CO2a = (cap/cMax)*CO2;



    //Solve for vapor transport
    fvScalarMatrix CO2Eqn =
    (
        (1 + (cap/cMax))*fvm::ddt(CO2)
      + fvm::div(phi,CO2)
      ==
        fvm::laplacian(((1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))))*D0/(3*tau)),CO2)
    );

    CO2Eqn.relax();

    fvConstraints().constrain(CO2Eqn);

    CO2Eqn.solve("CO2");

    fvConstraints().constrain(CO2);



    CO2a = (cap/cMax)*CO2;



    volScalarField diffCO2 = CO2 - CO2.oldTime();
    dCO2Avg = max(small,gAverage(diffCO2));
    dCO2Num = max(small,max(gMax(diffCO2),-gMin(diffCO2)));

    Info<< "CO2 Difference mean: " << dCO2Avg << " max: " << dCO2Num << endl;
}


// ************************************************************************* //
