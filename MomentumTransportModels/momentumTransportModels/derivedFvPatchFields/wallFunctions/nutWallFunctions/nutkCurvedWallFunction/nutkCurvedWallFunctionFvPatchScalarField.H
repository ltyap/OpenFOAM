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
    Foam::nutkCurvedWallFunctionFvPatchScalarField

Description
    Turbulent viscosity wall-function boundary condition
    for high Reynolds number flows based on near-wall turbulence kinetic energy.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutkCurvedWallFunction;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutkCurvedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkCurvedWallFunctionFvPatchScalarField_H
#define nutkCurvedWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nutkCurvedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkCurvedWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    // Private Data

        //- Curvature radius
        scalar Ri_;

        //- First constant for log-law  region
        scalar beta_;

        //- Second constant for log-law region
        scalar D_;
    protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutkCurvedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkCurvedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkCurvedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutkCurvedWallFunctionFvPatchScalarField
        //  onto a new patch
        nutkCurvedWallFunctionFvPatchScalarField
        (
            const nutkCurvedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nutkCurvedWallFunctionFvPatchScalarField
        (
            const nutkCurvedWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nutkCurvedWallFunctionFvPatchScalarField
        (
            const nutkCurvedWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutkCurvedWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions
            
        //- Return Ri
        scalar Ri() const
        {
            return Ri_;
        }

        //- Return beta
        scalar beta() const
        {
            return beta_;
        }

        //- Return D
        scalar D() const
        {
            return D_;
        }

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
