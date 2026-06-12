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
    Foam::nutUCurvedWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on velocity.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutUCurvedWallFunction;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutUCurvedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUCurvedWallFunctionFvPatchScalarField_H
#define nutUCurvedWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nutUCurvedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUCurvedWallFunctionFvPatchScalarField
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

        //- Calculate yPLus
        virtual tmp<scalarField> yPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutUCurvedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUCurvedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUCurvedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUCurvedWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUCurvedWallFunctionFvPatchScalarField
        (
            const nutUCurvedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nutUCurvedWallFunctionFvPatchScalarField
        (
            const nutUCurvedWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nutUCurvedWallFunctionFvPatchScalarField
        (
            const nutUCurvedWallFunctionFvPatchScalarField&,
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
                new nutUCurvedWallFunctionFvPatchScalarField(*this, iF)
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

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
