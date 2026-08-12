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

#include "fKmodel.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fKmodel, 0);
    defineRunTimeSelectionTable(fKmodel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fKmodel::fKmodel
(
    const word& name,
    const momentumTransportModel& turbulence
)
:
    momentumTransportModel_(turbulence)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fKmodel> Foam::fKmodel::New
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict
)
{
    const word fKmodelType(dict.lookup("fKmodel"));

    Info<< "Selecting fK model type " << fKmodelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fKmodelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fK model type "
            << fKmodelType << nl << nl
            << "Valid fK model types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fKmodel>(cstrIter()(name, turbulence, dict));
}


Foam::autoPtr<Foam::fKmodel> Foam::fKmodel::New
(
    const word& name,
    const momentumTransportModel& turbulence,
    const dictionary& dict,
    const dictionaryConstructorTable& additionalConstructors
)
{
    const word fKmodelType(dict.lookup("fKmodel"));

    Info<< "Selecting fK model type " << fKmodelType << endl;

    // First on additional ones
    dictionaryConstructorTable::const_iterator cstrIter =
        additionalConstructors.find(fKmodelType);

    if (cstrIter != additionalConstructors.end())
    {
        return autoPtr<fKmodel>(cstrIter()(name, turbulence, dict));
    }
    else
    {
        dictionaryConstructorTable::const_iterator cstrIter =
            dictionaryConstructorTablePtr_->find(fKmodelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown fK model type "
                << fKmodelType << nl << nl
                << "Valid fK model types are :" << endl
                << additionalConstructors.sortedToc()
                << " and "
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
            return autoPtr<fKmodel>();
        }
        else
        {
            return autoPtr<fKmodel>(cstrIter()(name, turbulence, dict));
        }
    }
}


// ************************************************************************* //
