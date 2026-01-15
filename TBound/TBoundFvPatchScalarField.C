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

#include "TBoundFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mappedPatchBase.H"

#include "fvcGrad.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::TBoundFvPatchScalarField::t() const
{
    return db().time().userTimeValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TBoundFvPatchScalarField::
TBoundFvPatchScalarField
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


Foam::TBoundFvPatchScalarField::
TBoundFvPatchScalarField
(
    const TBoundFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    dChanStrength_(ptf.dChanStrength_),
    boundSide_(ptf.boundSide_)
{}


Foam::TBoundFvPatchScalarField::
TBoundFvPatchScalarField
(
    const TBoundFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    dChanStrength_(ptf.dChanStrength_),
    boundSide_(ptf.boundSide_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TBoundFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);
}


void Foam::TBoundFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);
}


void Foam::TBoundFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr = refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    const word chanVal = "Tc";
    const word wallVal = "T";

    scalarField invDist = mpp.fromNeighbour(patchNbr.deltaCoeffs());

    if (chanVal == boundSide_)
    {
        //Also is center value for derivative condition
        scalarField dVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(wallVal).patchInternalField());

        scalarField sVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(wallVal));

        scalarField pVal = (sVal - dVal)*invDist;

	dimensionedScalar kappa_("kappa",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappa"));
	dimensionedScalar Xs_("Xs",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xs"));
	dimensionedScalar kappaS_("kappaS",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappaS"));
	dimensionedScalar Xa_("Xa",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xa"));
	dimensionedScalar kappaA_("kappaA",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappaA"));
	dimensionedScalar kappaL_("kappaL",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappaL"));

	scalar kappa = kappa_.value();
	scalar Xs = Xs_.value();
	scalar kappaS = kappaS_.value();
	scalar Xa = Xa_.value();
	scalar kappaA = kappaA_.value();
	scalar kappaL = kappaL_.value();

	scalarField Xl = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>("Xl"));

        mixedFvPatchScalarField::refGrad() = -(pVal)*((Xs*kappaS + Xa*kappaA + Xl*kappaL)/kappa);
        mixedFvPatchScalarField::refValue() = (sVal);

        valueFraction() = (dChanStrength_);
    }
    else if (wallVal == boundSide_)
    {
        //Also is center value for derivative condition
        scalarField dVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(chanVal).patchInternalField());

        scalarField sVal = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>(chanVal));

        scalarField pVal = (sVal - dVal)*invDist;

	dimensionedScalar kappa_("kappa",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappa"));
	dimensionedScalar Xs_("Xs",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xs"));
	dimensionedScalar kappaS_("kappaS",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappaS"));
	dimensionedScalar Xa_("Xa",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xa"));
	dimensionedScalar kappaA_("kappaA",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappaA"));
	dimensionedScalar kappaL_("kappaL",dimPower/(dimLength*dimTemperature),db().lookupObject<IOdictionary>("physicalProperties").lookup("kappaL"));

	scalar kappa = kappa_.value();
	scalar Xs = Xs_.value();
	scalar kappaS = kappaS_.value();
	scalar Xa = Xa_.value();
	scalar kappaA = kappaA_.value();
	scalar kappaL = kappaL_.value();

	scalarField Xl = patch().lookupPatchField<volScalarField, scalar>("Xl");

        mixedFvPatchScalarField::refGrad() = -(pVal)*(kappa/(Xs*kappaS + Xa*kappaA + Xl*kappaL));
        mixedFvPatchScalarField::refValue() = (sVal);

        valueFraction() = (1 - dChanStrength_);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::TBoundFvPatchScalarField::write
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
        TBoundFvPatchScalarField
    );
}

// ************************************************************************* //
