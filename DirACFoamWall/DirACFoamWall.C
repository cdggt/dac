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
#include "localEulerDdtScheme.H"
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(DirACFoamWall, 0);
    addToRunTimeSelectionTable(solver, DirACFoamWall, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Used to set maximum time step based on advection
void Foam::solvers::DirACFoamWall::correctCoNum()
{
    fluidSolver::correctCoNum(phi);
}

// Used to set maximum time step based on diffusion
void Foam::solvers::DirACFoamWall::correctDiNum()
{
    const volScalarField::Internal DiNumvf
    (
        fvc::surfaceSum
        (
            mesh.magSf()
	   *kappaS_
           *mesh.surfaceInterpolation::deltaCoeffs()
        )()()
       /(mesh.V()*rhoS_*CpS_)
       *runTime.deltaT()
    );

    const scalar meanDiNum = gAverage(DiNumvf);
    const scalar maxDiNum = gMax(DiNumvf);

    Info<< "Temp Diffusion Number mean: " << meanDiNum
        << " max: " << maxDiNum << endl;

    DiNum = maxDiNum;

    const volScalarField::Internal DiNumvf2
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

    const scalar meanDiNum2 = gAverage(DiNumvf2);
    const scalar maxDiNum2 = gMax(DiNumvf2);

    Info<< "Vapor Diffusion Number mean: " << meanDiNum2
        << " max: " << maxDiNum2 << endl;

    DiNum = max(DiNum,maxDiNum2);
}


void Foam::solvers::DirACFoamWall::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//To add your own variable to the solver, copy the code blocks for one of the
//folowing that is closest to what you want to add. If it includes "IOobject",
//then it is a variable tha4t can be solved for. If not, it is a constant you
//provide in the physicalProperties file. "j" is the outlier, a volume field
//that is not initialized but is calculated during runtime and output 
Foam::solvers::DirACFoamWall::DirACFoamWall(fvMesh& mesh)
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

    T_
    (
        IOobject
	(
	    "T",
	    runTime.name(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    ),

    Xv_
    (
        IOobject
	(
	    "Xv",
	    runTime.name(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    ),

    Xl_
    (
        IOobject
	(
	    "Xl",
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

    CpS_
    (
      "CpS",
      dimEnergy/(dimMass*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("CpS")
    ),

    rhoS_
    (
      "rhoS",
      dimDensity,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("rhoS")
     ),

    kappaS_
    (
      "kappaS",
      dimPower/(dimLength*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("kappaS")
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

    CpA_
    (
      "CpA",
      dimEnergy/(dimMass*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("CpA")
    ),

    rhoA_
    (
      "rhoA",
      dimDensity,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("rhoA")
     ),

    kappaA_
    (
      "kappaA",
      dimPower/(dimLength*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("kappaA")
    ),

    CpL_
    (
      "CpL",
      dimEnergy/(dimMass*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("CpL")
    ),

    rhoL_
    (
      "rhoL",
      dimDensity,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("rhoL")
    ),

    kappaL_
    (
      "kappaL",
      dimPower/(dimLength*dimTemperature),
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("kappaL")
    ),

    rhoV_
    (
      "rhoV",
      dimDensity,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("rhoV")
    ),

    j_
    (
        IOobject
        (
            "j",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Xv_ * dimensionedScalar("jVal", dimless/dimTime, 0)
    ),

    Rw_
    (
      "Rw",
      dimGasConstant,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("Rw")
    ),

    Lw_
    (
      "Lw",
      dimEnergy/dimMass,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("Lw")
    ),

    fact_
    (
      "fact",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("fact")
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

    sig_
    (
      "sig",
      dimForce/dimLength,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("sig")
    ),

    AH_
    (
      "AH",
      dimEnergy,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("AH")
    ),

    eta_
    (
      "eta",
      dimDynamicViscosity,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("eta")
    ),

    A_
    (
      "A",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("A")
    ),

    B_
    (
      "B",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("B")
    ),

    C_
    (
      "C",
      dimless,
      db().lookupObject<IOdictionary>
      (
          "physicalProperties"
      ).lookup("C")
    ),

    DiNum(0),
    CoNum(0),

    MRF(mesh),

    //Repeat every variable after here for the public reference

    p(p_),
    U(U_),
    phi(phi_),
    T(T_),
    Xv(Xv_),
    Xl(Xl_),
    Xs(Xs_),
    CpS(CpS_),
    rhoS(rhoS_),
    kappaS(kappaS_),
    Xa(Xa_),
    CpA(CpA_),
    rhoA(rhoA_),
    kappaA(kappaA_),
    CpL(CpL_),
    rhoL(rhoL_),
    kappaL(kappaL_),
    rhoV(rhoV_),
    j(j_),
    Rw(Rw_),
    Lw(Lw_),
    fact(fact_),
    D0(D0_),
    tau(tau_),
    lam(lam_),
    r0(r0_),
    sig(sig_),
    AH(AH_),
    eta(eta_),
    A(A_),
    B(B_),
    C(C_)

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

        dXvAvg = maxXvDelta;
        dXvNum = maxXvDelta;

        dXlAvg = maxXlDelta;
        dXlNum = maxXlDelta;

        dTAvg = maxTDelta;
        dTNum = maxTDelta;
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

Foam::solvers::DirACFoamWall::~DirACFoamWall()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::DirACFoamWall::readControls()
{
    maxDi = runTime.controlDict().lookupOrDefault<scalar>("maxDi", 1.0);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    maxXvDelta = runTime.controlDict().lookupOrDefault<scalar>("maxXvDelta", 1.0);

    maxXlDelta = runTime.controlDict().lookupOrDefault<scalar>("maxXlDelta", 1.0);

    maxTDelta = runTime.controlDict().lookupOrDefault<scalar>("maxTDelta", 1.0);

    maxDeltaT_ =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", vGreat);

    correctPhi = pimple.dict().lookupOrDefault
     (
         "correctPhi",
         mesh.dynamic()
     );
}

Foam::scalar Foam::solvers::DirACFoamWall::maxDeltaT() const
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

    if (maxXvDelta > small)
    {
        deltaT = min(deltaT, maxXvDelta/dXvNum*runTime.deltaTValue());
    }

    if (maxXlDelta > small)
    {
        deltaT = min(deltaT, maxXlDelta/dXlNum*runTime.deltaTValue());
    }

    if (maxTDelta > small)
    {
        deltaT = min(deltaT, maxTDelta/dTNum*runTime.deltaTValue());
    }

    return deltaT;
}

void Foam::solvers::DirACFoamWall::preSolve()
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


void Foam::solvers::DirACFoamWall::prePredictor()
{
    if (pimple.predictTransport())
    {
        momentumTransport->predict();
    }
}


void Foam::solvers::DirACFoamWall::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUEqn.clear();
}


void Foam::solvers::DirACFoamWall::postCorrector()
{
    if (pimple.correctTransport())
    {
        viscosity->correct();
        momentumTransport->correct();
    }
}


void Foam::solvers::DirACFoamWall::postSolve()
{}


// ************************************************************************* //
