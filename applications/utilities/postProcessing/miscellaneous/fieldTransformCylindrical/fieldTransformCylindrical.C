/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    

Description


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include "IOobjectList.H"
#include "cylindricalTransform.H"
// #include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transform fields into cylindrical coordinate. Assumes that the axial direction is the x-axis."
    );

    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

    IOdictionary fieldTransformCylindricalDict
    (
        IOobject
        (
            "fieldTransformCylindricalDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    cylindricalTransform cylTrans(mesh, fieldTransformCylindricalDict);

    const wordList& fieldList(fieldTransformCylindricalDict.lookup("fields"));

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Transforming fields for time " << runTime.timeName() << endl;
        forAll(fieldList, fieldI)
        {
            const word fieldName = fieldList[fieldI];
            Info << "Transforming " << fieldName << endl;
            // cylTrans.transformFields<scalar>(fieldName, runTime.timeName());
            cylTrans.transformFields<vector>(fieldName, runTime.timeName());
            cylTrans.transformFields<symmTensor>(fieldName, runTime.timeName());
            cylTrans.transformFields<tensor>(fieldName, runTime.timeName());
            // cylTrans.transformFields<sphericalTensor>(fieldName, runTime.timeName());
        }

    }
    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //