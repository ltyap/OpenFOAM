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
    Foam::fKmodels::GirimajiFk

Description
    Calculation of fK via Girimaji's estimate. Details are found in:
    \verbatim
        Branislav Basara, Zoran Pavlovic, Sharath Girimaji,
        A new approach for the calculation of the cut-off resolution parameter 
        in bridging methods for turbulent flow simulation,
        International Journal of Heat and Fluid Flow,
        Volume 74, 2018.
    \endverbatim

SourceFiles
    GirimajiFk.C

\*---------------------------------------------------------------------------*/

#ifndef GirimajiFk_H
#define GirimajiFk_H

#include "fKmodel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fKmodels
{

/*---------------------------------------------------------------------------*\
                      Class GirimajiFk Declaration
\*---------------------------------------------------------------------------*/

class GirimajiFk
:
    public fKmodel
{
    // Private Data

        scalar Cpans_;

    // Private Member Functions


public:

    //- Runtime type information
    TypeName("GirimajiFk");


    // Constructors

        //- Construct from name, momentumTransportModel and dictionary
        GirimajiFk
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        GirimajiFk(const GirimajiFk&) = delete;


    //- Destructor
    virtual ~GirimajiFk()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);


        // Calculate the fK values
        tmp<volScalarField> calcFk(const volScalarField&);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const GirimajiFk&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fKmodels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
