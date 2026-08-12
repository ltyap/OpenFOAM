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

#include "ChaouatSchiestel09Fk.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fKmodels
{
    defineTypeNameAndDebug(ChaouatSchiestel09Fk, 0);
    addToRunTimeSelectionTable(fKmodel, ChaouatSchiestel09Fk, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fKmodels::ChaouatSchiestel09Fk::ChaouatSchiestel09Fk
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict
)
:
    fKmodel(name, turbulence),
    beta_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "beta_",
            0.0495
        )
    ),
    eta_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "eta_",
            2.0/9.0
        )
    )
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fKmodels::ChaouatSchiestel09Fk::read(const dictionary& dict)
{
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<scalar>
    (
        "beta",
        beta_
    );
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<scalar>
    (
        "eta",
        eta_
    );

}

Foam::tmp<Foam::volScalarField> 
Foam::fKmodels::ChaouatSchiestel09Fk::calcFk(const volScalarField& delta)
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

    volScalarField piTimesLtByDeltaCubed = pow3(Foam::constant::mathematical::pi*lt/delta);

    tmp<volScalarField> tfK = 1/pow(1 + beta_*piTimesLtByDeltaCubed, eta_);

    return tfK;
}


// ************************************************************************* //
