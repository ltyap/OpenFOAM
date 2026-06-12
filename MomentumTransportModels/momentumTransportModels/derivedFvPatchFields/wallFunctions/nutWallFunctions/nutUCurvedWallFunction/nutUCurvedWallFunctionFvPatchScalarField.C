/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "nutUCurvedWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutUCurvedWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();
    //Info << "This is nutUCurvedWallFunction, function nut()" << endl;
    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& y = turbModel.y()[patchi];
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const scalarField yPlus(this->yPlus(magUp));

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    forAll(yPlus, facei)
    {
        scalar ratio = y[facei]/Ri_;
        //Info << y[facei] << endl;
        scalar g = (1/beta_)*log(y[facei]/(y[facei]+Ri_))+D_;
        if (yPlus[facei] > yPlusLam_)
        {
            nutw[facei] =
                nuw[facei]*(yPlus[facei]/g - 1);
        }
        else
        {
            nutw[facei] = nuw[facei]*(y[facei]/(Ri_*log(ratio + 1.0)) - 1.0);
            //Info << nutw[facei] << endl;
        }
    }
    // Info << "This is nutUCurvedWallFunction, function nut() end" << endl;

    return tnutw;
}


tmp<scalarField> nutUCurvedWallFunctionFvPatchScalarField::yPlus
(
    const scalarField& magUp
) const
{
    const label patchi = patch().index();
    //Info << "This is nutUCurvedWallFunction, function yPlus()" << endl;

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    //const vectorField& cf = patch().Cf();

    //Info << sqrt(sqr(cf.component(0))+ sqr(cf.component(1))) <<endl;
    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        const scalar Re = magUp[facei]*y[facei]/nuw[facei];
        const scalar ryPlusLam = 1/yPlusLam_;
        const scalar f = 1/(Ri_*log((y[facei]/Ri_)+1));
        const scalar g = (1/beta_)*log(y[facei]/(y[facei]+Ri_))+D_;

        int iter = 0;
        scalar yp = yPlusLam_;
        scalar yPlusLast = yp;

        do
        {
            yPlusLast = yp;
            if (yp > yPlusLam_)
            {
                yp = Re/g;
            }
            else
            { 
                yp = sqrt(y[facei]*Re*f);
            }
        } while(mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 20);

        yPlus[facei] = yp;
    }
    // Info << "This is nutUCurvedWallFunction, function yPlus() end" << endl;

    return tyPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUCurvedWallFunctionFvPatchScalarField::nutUCurvedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    Ri_(0.01),
    beta_(0.41),
    D_(5.5)
{}


nutUCurvedWallFunctionFvPatchScalarField::nutUCurvedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    Ri_(dict.lookup<scalar>("Ri")),
    beta_(dict.lookup<scalar>("beta")),
    D_(dict.lookup<scalar>("D"))
{}


nutUCurvedWallFunctionFvPatchScalarField::nutUCurvedWallFunctionFvPatchScalarField
(
    const nutUCurvedWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Ri_(ptf.Ri_),
    beta_(ptf.beta_),
    D_(ptf.D_)
{}


nutUCurvedWallFunctionFvPatchScalarField::nutUCurvedWallFunctionFvPatchScalarField
(
    const nutUCurvedWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf, iF),
    Ri_(sawfpsf.Cmu_),
    beta_(sawfpsf.beta_),
    D_(sawfpsf.D_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUCurvedWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    return yPlus(magUp);
}



void nutUCurvedWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "Ri", Ri_);
    writeEntry(os, "beta", beta_);
    writeEntry(os, "D", D_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUCurvedWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
