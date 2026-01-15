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

Class
    Foam::solvers::DirACFoamWallReduce

Description
    Solver module for steady or transient turbulent flow of incompressible
    isothermal fluids with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, constraining or limiting
    the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    DirACFoamWallReduce.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

My Notes
    A lot of this is copied over, so there may be some notes from the more
    general cases, especially above. I've done my best throughout the code
    to add comments that would be useful for someone familiar with C enough
    to read this code but not to work with it.

\*---------------------------------------------------------------------------*/

#ifndef DirACFoamWallReduce_H
#define DirACFoamWallReduce_H

#include "fluidSolver.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "pressureReference.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class DirACFoamWallReduce Declaration
\*---------------------------------------------------------------------------*/

class DirACFoamWallReduce
:
    public fluidSolver
{

protected:

    // Pressure

        //- Pressure field
        volScalarField p_;

        //- Pressure reference
        Foam::pressureReference pressureReference;


    // Kinematic properties

        //- Velocity field
        volVectorField U_;

        //- Volumetric-flux field
        surfaceScalarField phi_;


    // Momentum transport

        //- Kinematic viscosity model
        autoPtr<viscosityModel> viscosity;

        //- Pointer to the momentum transport model
        autoPtr<incompressible::momentumTransportModel> momentumTransport;

    // Temperature variable

        //- Temperature
        volScalarField T_;

    // Solid Material Information

        //- solid pore fraction
        dimensionedScalar Xs_;

        //- Heat Capacity solid
        dimensionedScalar CpS_;

        //- Density solid
        dimensionedScalar rhoS_;

        //- Thermal Conductivity solid (kappa to avoid conflict with other k terms)
        dimensionedScalar kappaS_;

    // Adsorbent Material Information

        //- adsorbent pore fraction
        dimensionedScalar Xa_;

        //- Heat Capacity adsorbent
        dimensionedScalar CpA_;

        //- Density adsorbent
        dimensionedScalar rhoA_;

        //- Thermal Conductivity adsorbent (kappa to avoid conflict with other k terms)
        dimensionedScalar kappaA_;

    // Liquid Water Heat Information

        //- Liquid concentration
        volScalarField Xl_;

        //- Heat Capacity of water
        dimensionedScalar CpL_;

        //- density of water
        dimensionedScalar rhoL_;

        //- Thermal Conductivity water (kappa to avoid conflict with other k terms)
        dimensionedScalar kappaL_;

    // Vapor Water Heat Information

        //- Vapor concentration
        volScalarField Xv_;

        //- density of vapor
        dimensionedScalar rhoV_;

    // Water Phase Change Information

        //- j for visualization in the output
        volScalarField j_;

        //- Specific Gas constant of water
        dimensionedScalar Rw_;

        //- Latent Heat of water
        dimensionedScalar Lw_;

        //- Artificial prefactor
        dimensionedScalar fact_;

    // Transport Variables

        //- Vapor Diffusivity
        dimensionedScalar D0_;

        //- Toruosity
        dimensionedScalar tau_;

        //- mean free path of vapor
        dimensionedScalar lam_;

        //- dry pore size
        dimensionedScalar r0_;

        //- surface tension
        dimensionedScalar sig_;

        //- Hamaker constant
        dimensionedScalar AH_;

        //- dynamic viscosity
        dimensionedScalar eta_;

    // Saturation Temperature

        //- A value
        dimensionedScalar A_;

        //- B value
        dimensionedScalar B_;

        //- C value
        dimensionedScalar C_;

    //Diffusive Time Step Control

        //- Diffusion number
        scalar DiNum;

        //- Courant number
        scalar CoNum;

        //- Max Xl change magnitude
        scalar dXlNum;

        //- Average Xl change
        scalar dXlAvg;

        //- Max T change magnitude
        scalar dTNum;

        //- Average T change
        scalar dTAvg;

        //- Max Diffusion number
        scalar maxDi;

        //- Max Courant number
        scalar maxCo;

        //- Max Xl change rate
        scalar maxXlDelta;

        //- Max T change rate
        scalar maxTDelta;

        //- max time step
        scalar maxDeltaT_;


    // Optional models

        //- MRF zone list
        IOMRFZoneList MRF;


    // Cached temporary fields

        //- Pointer to the surface momentum field
        //  used to recreate the flux after mesh-change
        autoPtr<surfaceVectorField> Uf;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;


    // Protected Member Functions

        //- Set rDeltaT for LTS
        void setRDeltaT();

        //- Read controls (my function)
        virtual void readControls();

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Correct the cached Diffusion numbers (my function)
        void correctDiNum();

        //- Calculate and print the continuity errors
        void continuityErrors();

        //- Construct the pressure equation
        //  and correct the pressure and velocity
        virtual void correctPressure();


public:

    // Public Data (The existing solvers liked using these for solving
    // equatons, so I followed suite here)

        //- Reference to the pressure field
        const volScalarField& p;

        //- Reference to the velocity field
        const volVectorField& U;

        //- Reference to the volumetric-flux field
        const surfaceScalarField& phi;

        //- Reference to density
        const volScalarField& T;

        //- Reference to solid pore fraction
        const dimensionedScalar& Xs;

        //- Reference to heat capacity
        const dimensionedScalar& CpS;

        //- Reference to density
        const dimensionedScalar& rhoS;

        //- Reference to thermal conductivity
        const dimensionedScalar& kappaS;

        //- Reference to adsorbent pore fraction
        const dimensionedScalar& Xa;

        //- Reference to heat capacity
        const dimensionedScalar& CpA;

        //- Reference to density
        const dimensionedScalar& rhoA;

        //- Reference to thermal conductivity
        const dimensionedScalar& kappaA;

        //- Reference to liquid concentration
        const volScalarField& Xl;

        //- Reference to heat capacity
        const dimensionedScalar& CpL;

        //- Reference to density of water
        const dimensionedScalar& rhoL;

        //- Reference to thermal conductivity
        const dimensionedScalar& kappaL;

        //- Reference to vapor concentration
        const volScalarField& Xv;

        //- Reference to density of vapor
        const dimensionedScalar& rhoV;

        //- Reference to j_
        const volScalarField& j;

        //- Reference to Specific Gas constant of water
        const dimensionedScalar& Rw;

        //- Reference to Latent Heat of water
        const dimensionedScalar& Lw;

        //- Reference to fact_
        const dimensionedScalar& fact;

        //- Reference to vapor diffusivity
        const dimensionedScalar& D0;

        //- Reference to tortuosity
        const dimensionedScalar& tau;

        //- Reference to mean free path
        const dimensionedScalar& lam;

        //- Reference to dry pore size
        const dimensionedScalar& r0;

        //- Reference to surface tension
        const dimensionedScalar& sig;

        //- Reference to Hamaker Constant;
        const dimensionedScalar& AH;

        //- Reference to dynamic viscosity;
        const dimensionedScalar& eta;

        //- Reference to A
        const dimensionedScalar& A;

        //- Reference to B
        const dimensionedScalar& B;

        //- Reference to C
        const dimensionedScalar& C;


    //- Runtime type information
    TypeName("DirACFoamWallReduce");


    // Constructors

        //- Construct from region mesh
        DirACFoamWallReduce(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        DirACFoamWallReduce(const DirACFoamWallReduce&) = delete;


    //- Destructor
    virtual ~DirACFoamWallReduce();


    // Member Functions

        //- Return maximum time step for stale solution (my solver)
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature (my version solves as temp)
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DirACFoamWallReduce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
