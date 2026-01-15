/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::monBoundSimpleFvPatchScalarField

Description
    This boundary condition provides a monBoundSimple condition,
    calculated as some ratio between value continuity and a flux continuity
    modified by porosity and tortuosity

Usage
    \table
        Property      | Description                           | Req'd?
        dChanStrength | strength of dirichlet on channel side | yes
        tau           | tortuosity of wall                    | yes
        phi           | surface porosity of wall              | yes
        boundSide     | Specify side of boundary              | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            monBoundSimple;
        dChanStrength   0.50;
        tau             4;
        phi             0.33333333;
        boundSide       O2;
    }
    \endverbatim

SourceFiles
    monBoundSimpleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef BC_H
#define BC_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class monBoundSimpleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class monBoundSimpleFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Single valued scalar quantity, e.g. a coefficient
        scalar dChanStrength_;

        //- String for what side of the condition we are on
        const word boundSide_;


    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("monBoundSimple");


    // Constructors

        //- Construct from patch, internal field and dictionary
        monBoundSimpleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField
        //  onto a new patch
        monBoundSimpleFvPatchScalarField
        (
            const monBoundSimpleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        monBoundSimpleFvPatchScalarField
        (
            const monBoundSimpleFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        monBoundSimpleFvPatchScalarField
        (
            const monBoundSimpleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new monBoundSimpleFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchScalarField&,
                const fvPatchFieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
