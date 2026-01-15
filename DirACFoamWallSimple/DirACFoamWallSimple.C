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
#include "localEulerDdtScheme.H"
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(DirACFoamWallSimple, 0);
    addToRunTimeSelectionTable(solver, DirACFoamWallSimple, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Used to set maximum time step based on advection
void Foam::solvers::DirACFoamWallSimple::correctCoNum()
{
    fluidSolver::correctCoNum(phi);
}

// Used to set maximum time step based on diffusion
void Foam::solvers::DirACFoamWallSimple::correctDiNum()
{
    const volScalarField::Internal DiNumvf
    (
        fvc::surfaceSum
        (
            mesh.magSf()
           *D0_
           *mesh.surfaceInterpolation::deltaCoeffs()
        )()()
       /(mesh.V())
       *runTime.deltaT()
    );

    const scalar meanDiNum = gAverage(DiNumvf);
    const scalar maxDiNum = gMax(DiNumvf);

    Info<< "Vapor Diffusion Number mean: " << meanDiNum
        << " max: " << maxDiNum << endl;

    DiNum = maxDiNum;
}


void Foam::solvers::DirACFoamWallSimple::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//To add your own variable to the solver, copy the code blocks for one of the
//folowing that is closest to what you want to add. If it includes "IOobject",
//then it is a variable tha4t can be solved for. If not, it is a constant you
//provide in the physicalProperties file. "j" is the outlier, a volume field
//that is not initialized but is calculated during runtime and output 
Foam::solvers::DirACFoamWallSimple::DirACFoamWallSimple(fvMesh& mesh)
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

    O2_
    (
        IOobject
	(
	    "O2",
	    runTime.name(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    ),

    Xs_
    (
      "Xs",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("Xs")
    ),

    Xa_
    (
      "Xa",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("Xa")
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

    tau_
    (
      "tau",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("tau")
    ),

    lam_
    (
      "lam",
      dimLength,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("lam")
    ),

    r0_
    (
      "r0",
      dimLength,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("r0")
    ),

    Xl_
    (
      "Xl",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("Xl")
    ),

    DiNum(0),
    CoNum(0),

    MRF(mesh),

    //Repeat every variable after here for the public reference

    p(p_),
    U(U_),
    phi(phi_),
    O2(O2_),
    Xs(Xs_),
    Xa(Xa_),
    D0(D0_),
    tau(tau_),
    lam(lam_),
    r0(r0_),
    Xl(Xl_)

//The only part you should touch after this part is the first if statement
//if you want to mess with the time step controls. Do not touch otherwise.
//(I did not touch anything else)
{
    readControls();

    mesh.schemes().setFluxRequired(p.name());

    momentumTransport->validate();

    if (transient())
    {
        correctCoNum();
        correctDiNum();

        dO2Avg = maxO2Delta;
        dO2Num = maxO2Delta;
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

//Don't touch anything after this part except for "readControls" if you again
//need to mess with time step control or if you know what you are doing and
//need to replace one of the functons.

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::DirACFoamWallSimple::~DirACFoamWallSimple()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamWallSimple::readControls()
{
    maxDi = runTime.controlDict().lookupOrDefault<scalar>("maxDi", 1.0);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxO2Delta = runTime.controlDict().lookupOrDefault<scalar>("maxO2Delta", 1.0);

    maxDeltaT_ =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", vGreat);

    correctPhi = pimple.dict().lookupOrDefault
     (
         "correctPhi",
         mesh.dynamic()
     );
}

Foam::scalar Foam::solvers::DirACFoamWallSimple::maxDeltaT() const
{

    scalar deltaT = min(fvModels().maxDeltaT(), maxDeltaT_);

    if (CoNum > small)
    {
        deltaT = min(deltaT, maxCo/CoNum*runTime.deltaTValue());
    }

    if (DiNum > small)
    {
        deltaT = min(deltaT, maxDi/DiNum*runTime.deltaTValue());
    }

    if (maxO2Delta > small)
    {
        deltaT = min(deltaT, maxO2Delta/dO2Num*runTime.deltaTValue());
    }

    return deltaT;
}

void Foam::solvers::DirACFoamWallSimple::preSolve()
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
        correctDiNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::DirACFoamWallSimple::prePredictor()
{
    if (pimple.predictTransport())
    {
        momentumTransport->predict();
    }
}


void Foam::solvers::DirACFoamWallSimple::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUEqn.clear();
}


void Foam::solvers::DirACFoamWallSimple::postCorrector()
{
    if (pimple.correctTransport())
    {
        viscosity->correct();
        momentumTransport->correct();
    }
}


void Foam::solvers::DirACFoamWallSimple::postSolve()
{}


// ************************************************************************* //
