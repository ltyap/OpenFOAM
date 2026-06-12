/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |               2020 Timofey Mukha
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

Application
    postChannel

Description
    Post-processes data from channel flow calculations.

    Assuming that the mesh is periodic in the x and z directions, collapse
    fields to a line and print them to postProcesing/collapsedFields.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "channelIndex.H"
#include "makeGraph.H"

#include "OSspecific.H"
#include "IOobjectList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
Foam::wordList comptNames()
{
    return wordList();
}


template<>
Foam::wordList comptNames<vector>()
{
    return {"_X", "_Y", "_Z"};
}


template<>
Foam::wordList comptNames<symmTensor>()
{
    return {"_XX", "_XY", "_XZ", "_YY", "_YZ", "_ZZ"};
}


template<>
Foam::wordList comptNames<tensor>()
{
    return {"_XX", "_XY", "_XZ", "_YX", "_YY", "_YZ", "_ZX", "_ZY", "_ZZ"};
}

template<>
Foam::wordList comptNames<sphericalTensor>()
{
    return {"_II"};
}


template<class T>
void writeToFile
(
    const scalarField& y,
    const Field<T>& values,
    const word name,
    const fileName path,
    const word format,
    const wordList comptNames)
{
    for (label i=0; i<comptNames.size(); ++i)
    {
        makeGraph(y, values.component(i), name+comptNames[i], path, format);
    }
}

template<>
void writeToFile
(
    const scalarField& y,
    const Field<scalar>& values,
    const word name,
    const fileName path,
    const word format,
    const wordList comptNames)
{
    makeGraph(y, values, name, path, format);
}

template<class T>
void collapse
(
    const word fieldName,
    const word timeName,
    const fvMesh & mesh,
    const channelIndex & channelIndexing,
    const word format
)
{
    bool debug = false;
    using FieldType = GeometricField<T, fvPatchField, volMesh>;

    IOobject fieldHeader
    (
        fieldName,
        timeName,
        mesh,
        IOobject::MUST_READ
    );

    if (debug)
    {
        Info << "Initial type name of fieldHeader: "
        <<fieldHeader.headerClassName() << endl;
    }
    fileName path(fieldHeader.rootPath()/fieldHeader.caseName()/
    "postProcessing"/"collapsedFields"/fieldHeader.instance());

    // There is no option to suppress the warning message
    // in the function typeHeaderOk for OpenFOAM v9. 
    // However, we also need to run it to update the 
    // variable headerClassName_ in the class IOobject. 
    // Hence we set the flag to false and then 
    // do the check ourselves.
    fieldHeader.typeHeaderOk<FieldType>(false);
    word suppliedTypeName = FieldType::typeName;
    word headerClassName = fieldHeader.headerClassName();
    if (debug)
    {
        Info << "Type name as read from input file: " << headerClassName << endl;
        Info << "Type name as supplied by template: "<< suppliedTypeName << endl;
    }

    if (suppliedTypeName == headerClassName)
    {

        mkDir(path);
        FieldType field
        (
            fieldHeader,
            mesh
        );

        Pair<T> boundaryValues =
            channelIndexing.collapseBoundary<T>
            (
                mesh.boundaryMesh(),
                field.boundaryField()
            );

        Field<T> internalValues =
            channelIndexing.collapse(field);

        Field<T> allValues(internalValues.size() + 2);
        allValues[0] = boundaryValues[0];
        allValues[allValues.size()-1] = boundaryValues[1];

        for (label i=0; i<internalValues.size(); ++i)
        {
            allValues[i+1] = internalValues[i];
        }

        tmp<scalarField> yTmp = channelIndexing.y(mesh.boundaryMesh());
        scalarField& y = yTmp.ref();

        const wordList cNames = comptNames<T>();
        writeToFile<T>(y, allValues, fieldName, path, format, cNames);
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Post-process data from channel flow calculations"
    );

    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"
#   include "readTransportProperties.H"

    const word& gFormat = runTime.graphFormat();

    // Setup channel indexing for averaging over channel down to a line

    IOdictionary annularDict
    (
        IOobject
        (
            "postAnnularDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    channelIndex channelInd(mesh, annularDict);

    const wordList& fieldList(annularDict.lookup("fields"));

 
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Collapsing fields for time " << runTime.timeName() << endl;

        forAll(fieldList, fieldI)
        {
            const word fieldName = fieldList[fieldI];
            Info << "Collapsing " << fieldName << endl;
            collapse<scalar>(fieldName,runTime.timeName(), mesh, channelInd, gFormat);
            collapse<vector>(fieldName,runTime.timeName(), mesh, channelInd, gFormat);
            collapse<sphericalTensor>(fieldName,runTime.timeName(), mesh, channelInd, gFormat);
            collapse<symmTensor>(fieldName,runTime.timeName(), mesh, channelInd, gFormat);
            collapse<tensor>(fieldName,runTime.timeName(), mesh, channelInd, gFormat);
        }
    }
/* */
    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
