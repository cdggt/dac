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

\*---------------------------------------------------------------------------*/

#include "monBoundFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::monBoundFvPatchScalarField::t() const
{
    return db().time().userTimeValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::monBoundFvPatchScalarField::
monBoundFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    dChanStrength_(dict.lookup<scalar>("dChanStrength")),
    boundSide_(dict.lookup<word>("boundSide"))
{
    refValue() = patchInternalField();
    refGrad() = Zero;
    valueFraction() = 0.0;

    fvPatchScalarField::operator=(refValue());
}


Foam::monBoundFvPatchScalarField::
monBoundFvPatchScalarField
(
    const monBoundFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    dChanStrength_(ptf.dChanStrength_),
    boundSide_(ptf.boundSide_)
{}


Foam::monBoundFvPatchScalarField::
monBoundFvPatchScalarField
(
    const monBoundFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    dChanStrength_(ptf.dChanStrength_),
    boundSide_(ptf.boundSide_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::monBoundFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);
}


void Foam::monBoundFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);
}


void Foam::monBoundFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr = refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    const word chanVal = "Xc";
    const word wallVal = "Xv";

    //scalarField invDist = mpp.fromNeighbour(patch().deltaCoeffs());

    scalarField invDist = mpp.fromNeighbour(patchNbr.deltaCoeffs());

    if (chanVal == boundSide_)
    {
        //Also is center value for derivative condition
        scalarField dVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(wallVal).patchInternalField());

        scalarField sVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(wallVal));

        scalarField pVal = (sVal - dVal)*invDist;

	dimensionedScalar Xs_("Xs",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xs"));
	dimensionedScalar Xa_("Xa",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xa"));
	dimensionedScalar lam_("lam",dimLength,db().lookupObject<IOdictionary>("physicalProperties").lookup("lam"));
	dimensionedScalar r0_("r0",dimLength,db().lookupObject<IOdictionary>("physicalProperties").lookup("r0"));
	dimensionedScalar tau_("tau",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("tau"));

	scalar Xs = Xs_.value();
	scalar Xa = Xa_.value();
	scalar lam = lam_.value();
	scalar r0 = r0_.value();
	scalar tau = tau_.value();

	scalarField Xl = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>("Xl"));

        mixedFvPatchScalarField::refGrad() = -((1-Xs-Xa-Xl)/(3*tau))*(1/(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa))))))*(pVal);
        mixedFvPatchScalarField::refValue() = (sVal);

        valueFraction() = (dChanStrength_);
    }
    else if (wallVal == boundSide_)
    {
        //Also is center value for derivative condition
        scalarField dVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(chanVal).patchInternalField());

        scalarField sVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(chanVal));

        scalarField pVal = (sVal - dVal)*invDist;

	dimensionedScalar Xs_("Xs",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xs"));
	dimensionedScalar Xa_("Xa",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xa"));
	dimensionedScalar lam_("lam",dimLength,db().lookupObject<IOdictionary>("physicalProperties").lookup("lam"));
	dimensionedScalar r0_("r0",dimLength,db().lookupObject<IOdictionary>("physicalProperties").lookup("r0"));
	dimensionedScalar tau_("tau",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("tau"));

	scalar Xs = Xs_.value();
	scalar Xa = Xa_.value();
	scalar lam = lam_.value();
	scalar r0 = r0_.value();
	scalar tau = tau_.value();

	scalarField Xl = patch().lookupPatchField<volScalarField, scalar>("Xl");

        mixedFvPatchScalarField::refGrad() = -((3*tau)/(1-Xs-Xa-Xl))*(1 + (lam/(r0*sqrt((1-Xs-Xa-Xl)/(1-Xs-Xa)))))*(pVal);
        mixedFvPatchScalarField::refValue() = (sVal);

        valueFraction() = (1 - dChanStrength_);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::monBoundFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "dChanStrength", dChanStrength_);
    writeEntry(os, "boundSide", boundSide_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        monBoundFvPatchScalarField
    );
}

// ************************************************************************* //
