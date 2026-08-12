/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "GirimajiFk.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fKmodels
{
    defineTypeNameAndDebug(GirimajiFk, 0);
    addToRunTimeSelectionTable(fKmodel, GirimajiFk, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fKmodels::GirimajiFk::GirimajiFk
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict
)
:
    fKmodel(name, turbulence),
    Cpans_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "Cpans",
            1/sqrt(0.09)
        )
    )
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fKmodels::GirimajiFk::read(const dictionary& dict)
{
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<scalar>
    (
        "Cpans",
        Cpans_
    );

}

Foam::tmp<Foam::volScalarField> 
Foam::fKmodels::GirimajiFk::calcFk(const volScalarField& delta)
{
    const fvMesh& mesh = momentumTransportModel_.mesh();

    const tmp<volScalarField> tk = momentumTransportModel_.k();
    const volScalarField& k = tk();
    
    // epsilon calculated as betaStar*omega*k for k-omega based models
    const tmp<volScalarField> teps = momentumTransportModel_.epsilon();
    const volScalarField& eps = teps();

    const volScalarField& kSSV = mesh.template lookupObject<volScalarField>("kSSV");


    const bool useOmega{mesh.template foundObject<volScalarField>("omega")};

    // Assume only either omega or epsilon is used and exists in the directory
    volScalarField lt = useOmega? sqrt(k + kSSV)/(eps/k) : pow(k + kSSV, 1.5)/eps;

    tmp<volScalarField> tfK = min
    (
        Cpans_*pow(delta/lt, 2.0/3.0),
        scalar(1.0)
    );

    return tfK;
    // return Foam::volScalarField::New
    // (
    //     "fK",

    // );
}


// ************************************************************************* //
