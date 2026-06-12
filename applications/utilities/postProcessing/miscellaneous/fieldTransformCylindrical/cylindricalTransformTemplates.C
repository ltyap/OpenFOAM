 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | Website:  https://openfoam.org
     \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
 
 #include "cylindricalTransform.H"
 #include "volFields.H"
 #include "surfaceFields.H"
 #include "transformGeometricField.H"
//  #include "transform.H"
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
template<class FieldType>
void Foam::cylindricalTransform::transformField
(
    const FieldType& field,
    const word timeName
)
{
    word transFieldName(transformFieldName(field.name()));
    Info << transFieldName << endl;

	FieldType transField
	(
		IOobject
        (   
            transFieldName,
            timeName,
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Foam::transform(RotationTensorFieldX(), field)
	);
    transField.write();
    // Info << transField << endl;
}


template<class Type>
void Foam::cylindricalTransform::transformFields
(
    const word fieldName,
    const word timeName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    // typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    IOobject fieldHeader
    (
        fieldName,
        timeName,
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    
    if
    (
        fieldHeader.typeHeaderOk<VolFieldType>(false)
    && fieldHeader.headerClassName() == VolFieldType::typeName
    )
    {
        // DebugInfo
        //     << type() << ": Field " << fieldName << " read from file"
        //     << endl;
        VolFieldType field
        (
            fieldHeader,
            mesh_
        );

        transformField<VolFieldType>(field, timeName);

    }
    // else if
    // (
    //     fieldHeader.typeHeaderOk<SurfaceFieldType>(false)
    // && fieldHeader.headerClassName() == SurfaceFieldType::typeName
    // )
    // {
    //     // DebugInfo
    //     //     << type() << ": Field " << fieldName << " read from file"
    //     //     << endl;

    //     transformField<SurfaceFieldType>
    //     (
    //         mesh_.lookupObject<SurfaceFieldType>(fieldName)
    //     );
    // }
    
}
