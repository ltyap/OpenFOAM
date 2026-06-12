/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::cylindricalTransform

Description


SourceFiles
    cylindricalTransform.C

\*---------------------------------------------------------------------------*/
#ifndef cylindricalTransform_H
#define cylindricalTransform_H

// #include "scalarField.H"
// #include "polyMesh.H"
#include "GeometricField.H"
#include "volMesh.H"
#include "fvPatchField.H"


namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class cylindricalTransform Declaration
\*---------------------------------------------------------------------------*/
class cylindricalTransform
{
    // Private data

        // axis of rotation
        const vector axis_;
        // origin
        const vector origin_;
        //Rotation matrix
        mutable autoPtr<volTensorField> RxPtr_;

        const fvMesh& mesh_;

    // Private Member Functions    
        tensor R(const vector& p) const;

        tmp<tensorField> RotX(const vectorField& cellCentre) const;

        void calcRotationTensorFieldX() const;


        //- No copy construct
        cylindricalTransform(const cylindricalTransform&) = delete;

        //- No copy assignment
        void operator=(const cylindricalTransform&) = delete;

public:

    //Constructors

        //- Construct from dictionary

        cylindricalTransform(const fvMesh&, const dictionary&);

    //Member functions

        //Return name of transformed field
        Foam::word transformFieldName(const word& fieldName) const;

        //- Transform the given field
        template<class FieldType>
        void transformField(const FieldType& field, const word timeName);
 
        //- Transform the given field if has the specified element type
        template<class Type>
        void transformFields(const word fieldName, const word timeName);

        //Access
        const volTensorField& RotationTensorFieldX()
        {
            return *RxPtr_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cylindricalTransformTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //