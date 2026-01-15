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
#include "localEulerDdtScheme.H"
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(DirACFoamChan, 0);
    addToRunTimeSelectionTable(solver, DirACFoamChan, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::DirACFoamChan::correctCoNum()
{
    fluidSolver::correctCoNum(phi);
}


void Foam::solvers::DirACFoamChan::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::DirACFoamChan::DirACFoamChan(fvMesh& mesh)
:
    fluidSolver(mesh),

    p_
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pressureReference(p_, pimple.dict()),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U_) & mesh.Sf()
    ),

    viscosity(viscosityModel::New(mesh)),

    momentumTransport
    (
        incompressible::momentumTransportModel::New
        (
            U_,
            phi_,
            viscosity
        )
    ),

    Tc_
    (
        IOobject
	(
	    "Tc",
	    runTime.name(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    ),

    Cp_
    (
      "Cp",
      dimEnergy/(dimMass*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("Cp")
    ),

    rho_
    (
      "rho",
      dimDensity,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("rho")
     ),

    kappa_
    (
      "kappa",
      dimPower/(dimLength*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("kappa")
    ),

    Xc_
    (
        IOobject
	(
	    "Xc",
	    runTime.name(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    ),

    D0_
    (
      "D0",
      dimArea/dimTime,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("D0")
    ),

    MRF(mesh),

    p(p_),
    U(U_),
    phi(phi_),
    Tc(Tc_),
    Cp(Cp_),
    rho(rho_),
    kappa(kappa_),
    Xc(Xc_),
    D0(D0_)
{
    mesh.schemes().setFluxRequired(p.name());

    momentumTransport->validate();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::DirACFoamChan::~DirACFoamChan()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamChan::preSolve()
{
    // Read the controls
    readControls();

    if ((mesh.dynamic() || MRF.size()) && !Uf.valid())
    {
        Info<< "Constructing face momentum Uf" << endl;

        // Ensure the U BCs are up-to-date before constructing Uf
        U_.correctBoundaryConditions();

        Uf = new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U)
        );
    }

    fvModels().preUpdateMesh();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::DirACFoamChan::prePredictor()
{
    if (pimple.predictTransport())
    {
        momentumTransport->predict();
    }
}


//void Foam::solvers::DirACFoamChan::thermophysicalPredictor()
//{}


void Foam::solvers::DirACFoamChan::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUEqn.clear();
}


void Foam::solvers::DirACFoamChan::postCorrector()
{
    if (pimple.correctTransport())
    {
        viscosity->correct();
        momentumTransport->correct();
    }
}


void Foam::solvers::DirACFoamChan::postSolve()
{}


// ************************************************************************* //
