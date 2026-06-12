/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "cylindricalTransform.H"
// #include "meshTools.H"
// #include "Time.H"
#include "axesRotation.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//Adapted from cylindrical.H
Foam::tensor Foam::cylindricalTransform::R(const vector& p) const
{
    vector dir{p - origin_};
    dir /= mag(dir) + vSmall;
    const vector r{dir - (dir & axis_)*axis_};
    
     if (mag(r) < small)
     {
        // If the point is on the axis choose any radial direction
        axesRotation ax{axis_, perpendicular(axis_)};
        return tensor(ax.e3(), ax.e1(), ax.e2());
     }
     else
     {
        axesRotation ax{axis_, dir};
        return tensor(ax.e3(), ax.e1(), ax.e2());
     }
}

Foam::tmp<Foam::tensorField> 
Foam::cylindricalTransform::RotX(const vectorField& cellCentre) const
{
    tmp<tensorField> rotTmp{new tensorField(cellCentre.size(), Foam::Zero)};
    tensorField& rot = rotTmp.ref();
    
    forAll(rot, cellI)
    {
        rot[cellI] = this->R(cellCentre[cellI]);

        /*
        const scalar z{cellCentre[cellI].component(2)};
        const scalar y{cellCentre[cellI].component(1)};
        scalar theta{Foam::atan2(z, y)};

        rot[cellI] = tensor
                        (
                            axis_.x(), axis_.y(), axis_.z(),
                            0, Foam::cos(theta), Foam::sin(theta),
                            0, -Foam::sin(theta), Foam::cos(theta)
                        );
        if (cellI == 1)
        {
            Info << "z: " << z << endl;
            Info <<"y: " <<  y << endl;
            Info << "theta: " << theta << endl;
            Info << "tensor: " << rot[cellI] << endl;
            Info << "---------" << endl;
        }
        // Info << "tensor: " << this->R(cellCentre[cellI]) << endl;
        */

    }
    return rotTmp;
}

void Foam::cylindricalTransform::calcRotationTensorFieldX
() const
{
    using FieldType = volTensorField;
    using BoundaryType = volTensorField::Boundary;

    const vectorField& cellcentres_ = mesh_.cellCentres();
    if (!RxPtr_.valid())
    {
        Info << "hi" << endl;

        tensorField rotations = RotX(cellcentres_);
        FieldType* volRotation = new FieldType
            (
                IOobject
                (
                    "volRotation",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedTensor(dimless, Zero)
            );

        RxPtr_.reset(volRotation);

        auto& vrotTensor = *RxPtr_;
        vrotTensor.primitiveFieldRef() = rotations;

        // Boundaries
        BoundaryType& bf = const_cast<BoundaryType&>(vrotTensor.boundaryField());
        forAll(bf, patchi)
        {
            bf[patchi] = RotX(bf[patchi].patch().patch().faceCentres());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::cylindricalTransform::cylindricalTransform
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    axis_(Foam::vector(1,0,0)),
    origin_(dict.lookupOrDefault<vector>("origin", Foam::vector(0,0,0))),
    RxPtr_(),
    mesh_(mesh)
{
    // const vectorField& cellCentres = mesh.cellCentres();

    calcRotationTensorFieldX();
    // Info << RotationTensorFieldX() << endl;
}

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

Foam::word Foam::cylindricalTransform::transformFieldName
(
    const word& fieldName
) const
{
    return fieldName + "Transformed";
}