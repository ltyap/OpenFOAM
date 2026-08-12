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

Class
    Foam::fKmodel

Description
    Abstract base class for PANS fK models

SourceFiles
    fKmodel.C

\*---------------------------------------------------------------------------*/

#ifndef fKmodel_H
#define fKmodel_H

#include "momentumTransportModel.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class fKmodel
{
protected:

    // Protected data

        const momentumTransportModel& momentumTransportModel_;


public:

    //- Runtime type information
    TypeName("fKmodel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fKmodel,
            dictionary,
            (
                const word& name,
                const momentumTransportModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, momentumTransportModel and dictionary
        fKmodel
        (
            const word& name,
            const momentumTransportModel& turbulence
        );

        //- Disallow default bitwise copy construction
        fKmodel(const fKmodel&) = delete;


    // Selectors

        //- Return a reference to the selected fK model
        static autoPtr<fKmodel> New
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary& dict
        );

        //- Return a reference to the selected fK model
        static autoPtr<fKmodel> New
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable&
        );


    //- Destructor
    virtual ~fKmodel()
    {}


    // Member Functions

        //- Return momentumTransportModel reference
        const momentumTransportModel& turbulence() const
        {
            return momentumTransportModel_;
        }

        //- Read the fKmodel dictionary
        virtual void read(const dictionary&) = 0;

        virtual tmp<volScalarField> calcFk(const volScalarField&) = 0;

    // Member Operators

        void operator=(const fKmodel&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
