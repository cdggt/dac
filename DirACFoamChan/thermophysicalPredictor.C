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

#include "DirACFoamChan.H"
#include "fvcGrad.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamChan::thermophysicalPredictor()
{
    volScalarField& Xc(Xc_);
    fvScalarMatrix XcEqn =
    (
        fvm::ddt(Xc)
      + fvm::div(phi,Xc)
      ==
        fvm::laplacian(D0,Xc)
    );

    XcEqn.relax();

    fvConstraints().constrain(XcEqn);

    XcEqn.solve("Xc");

    
    volScalarField& Tc(Tc_);

    fvScalarMatrix EEqn =
    (
        fvm::ddt(Cp*rho,Tc) + fvm::div(Cp*rho*phi,Tc)
     ==
        fvm::laplacian(kappa,Tc)
    );

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();
}


// ************************************************************************* //
