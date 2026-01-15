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

#include "waterBoundFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waterBoundFvPatchScalarField::t() const
{
    return db().time().userTimeValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waterBoundFvPatchScalarField::
waterBoundFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{


    //fixedValueFvPatchScalarField::evaluate();


    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        patchInternalField()
    );

}


Foam::waterBoundFvPatchScalarField::
waterBoundFvPatchScalarField
(
    const waterBoundFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::waterBoundFvPatchScalarField::
waterBoundFvPatchScalarField
(
    const waterBoundFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waterBoundFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::map(ptf, mapper);
}


void Foam::waterBoundFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);
}


void Foam::waterBoundFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr = refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    scalarField Xv = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>("Xc"));
    scalarField T  = mpp.fromNeighbour(patchNbr.lookupPatchField<volScalarField, scalar>("Tc"));

    dimensionedScalar AH_("AH",dimEnergy,db().lookupObject<IOdictionary>("physicalProperties").lookup("AH"));
    dimensionedScalar r0_("r0",dimLength,db().lookupObject<IOdictionary>("physicalProperties").lookup("r0"));
    dimensionedScalar sig_("sig",dimForce/dimLength,db().lookupObject<IOdictionary>("physicalProperties").lookup("sig"));
    dimensionedScalar rhoL_("rhoL",dimDensity,db().lookupObject<IOdictionary>("physicalProperties").lookup("rhoL"));
    dimensionedScalar Lw_("Lw",dimEnergy/dimMass,db().lookupObject<IOdictionary>("physicalProperties").lookup("Lw"));
    dimensionedScalar A_("A",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("A"));
    dimensionedScalar B_("B",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("B"));
    dimensionedScalar C_("C",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("C"));
    dimensionedScalar Xs_("Xs",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xs"));
    dimensionedScalar Xa_("Xa",dimless,db().lookupObject<IOdictionary>("physicalProperties").lookup("Xa"));

    scalar AH = AH_.value();
    scalar r0 = r0_.value();
    scalar sig = sig_.value();
    scalar rhoL = rhoL_.value();
    scalar Lw = Lw_.value();
    scalar A = A_.value();
    scalar B = B_.value();
    scalar C = C_.value();
    scalar Xs = Xs_.value();
    scalar Xa = Xa_.value();

    scalarField T1 = (AH)/(6*Foam::constant::mathematical::pi*pow(r0,3)) * (T/T);
    scalarField T2 = (sig)/(r0) * (Xv/Xv);
    scalarField T3 = rhoL * Lw * ((T/((B/(A-log10(101325*Xv)))-C))-1);

    scalarField Xi = Xv/Xv;

    scalarField A0 = T3/T3;
    scalarField B0 = -(3*T3+T2)/T3;
    scalarField C0 = 3*(T3 + T2)/T3;
    scalarField D0 = -(T3 + 3*T2 - T1)/T3;
    scalarField E0 = T2/T3;

    scalarField a = -((3*pow(B0,2))/(8*pow(A0,2))) + (C0/A0);
    scalarField b = ((pow(B0,3))/(8*pow(A0,3))) - (B0*C0/(2*pow(A0,2))) + (D0/A0);
    scalarField c = -((3*pow(B0,4))/(256*pow(A0,4))) + ((C0*pow(B0,2))/(16*pow(A0,3))) - ((B0*D0)/(4*pow(A0,2)))  + (E0/A0);

    scalarField P = -(pow(a,2)/12) - c;
    scalarField Q = -(pow(a,3)/108) + (a*c/3) - (pow(b,2)/8);
    scalarField R = -(Q/2) + sqrt((pow(Q,2)/4)+(pow(P,3)/27));

    scalarField U = cbrt(R);
    scalarField y = -(5*a/6) + U - (P/(3*U));
    scalarField W = sqrt(a + 2*y);

    Xi = -(B0/(4*A0)) + (-W + sqrt(-3*a - 2*y + (2*b/W)))/2;

    scalarField Xl = (1 - Xs - Xa) * (1 - pow(Xi,2));

    fixedValueFvPatchScalarField::operator==
    (
        Xl
    );


    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::waterBoundFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waterBoundFvPatchScalarField
    );
}

// ************************************************************************* //
