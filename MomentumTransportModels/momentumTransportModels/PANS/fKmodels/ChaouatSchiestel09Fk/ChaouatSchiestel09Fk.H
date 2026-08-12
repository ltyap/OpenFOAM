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
    Foam::fKmodels::ChaouatSchiestel09Fk

Description
    Calculation of fK via Chaouat and Schiestel's estimate from 2009. Details are found in:
    \verbatim
        Bruno Chaouat, Roland Schiestel,
        Progress in subgrid-scale transport modelling for continuous hybrid non-zonal RANS/LES simulations,
        International Journal of Heat and Fluid Flow,
        Volume 30, Issue 4, 2009.
    \endverbatim

SourceFiles
    ChaouatSchiestel09Fk.C

\*---------------------------------------------------------------------------*/

#ifndef ChaouatSchiestel09Fk_H
#define ChaouatSchiestel09Fk_H

#include "fKmodel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fKmodels
{

/*---------------------------------------------------------------------------*\
                      Class ChaouatSchiestel09Fk Declaration
\*---------------------------------------------------------------------------*/

class ChaouatSchiestel09Fk
:
    public fKmodel
{
    // Private Data

        scalar beta_;
        scalar eta_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("ChaouatSchiestel09Fk");


    // Constructors

        //- Construct from name, momentumTransportModel and dictionary
        ChaouatSchiestel09Fk
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        ChaouatSchiestel09Fk(const ChaouatSchiestel09Fk&) = delete;


    //- Destructor
    virtual ~ChaouatSchiestel09Fk()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);


        // Calculate the fK values
        tmp<volScalarField> calcFk(const volScalarField&);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ChaouatSchiestel09Fk&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fKmodels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
