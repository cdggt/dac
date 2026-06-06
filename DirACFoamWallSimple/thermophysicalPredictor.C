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

#include "DirACFoamWallSimple.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

#include "fvcSmooth.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamWallSimple::thermophysicalPredictor()
{

    //Make references to all volume field variables
    volScalarField& O2(O2_);



    //Solve for vapor transport
    fvScalarMatrix O2Eqn =
    (
        fvm::ddt(O2)
      + fvm::div(phi,O2)
      ==
        fvm::laplacian((1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))))*D0/(3*tau),O2)
    );

    O2Eqn.relax();

    fvConstraints().constrain(O2Eqn);

    O2Eqn.solve("O2");

    fvConstraints().constrain(O2);



    volScalarField diffO2 = fvc::laplacian((1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))))*D0/(3*tau),O2);

    dO2Avg = gAverage(diffO2)*runTime.deltaTValue();
    dO2Num = max(gMax(diffO2),-gMin(diffO2))*runTime.deltaTValue();

    Info<< "O2 Difference mean: " << dO2Avg << " max: " << dO2Num << endl;
}


// ************************************************************************* //
