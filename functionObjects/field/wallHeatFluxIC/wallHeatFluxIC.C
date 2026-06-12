/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "wallHeatFluxIC.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFluxIC, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFluxIC, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFluxIC::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall heat flux");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::wallHeatFluxIC::calcHeatFlux
(
    const volVectorField& gradT
)
{
    tmp<volScalarField> twallHeatFluxIC
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar(gradT.dimensions()*dimViscosity, Zero)
        )
    );
    if (!mesh_.foundObject<volScalarField>("alpha"))
    {
        FatalErrorInFunction
            << "Unable to find alpha field in the "
            << "database" << exit(FatalError);
    }

    const volScalarField& alpha = mesh_.lookupObjectRef<volScalarField>("alpha");

    volScalarField::Boundary& wallHeatFluxICBf =
        twallHeatFluxIC.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
       
        const scalarField& alphaw = alpha.boundaryField()[patchi];

        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        wallHeatFluxICBf[patchi] = alphaw*(-Sfp/magSfp) & gradT.boundaryField()[patchi];
    }

    return twallHeatFluxIC;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxIC::wallHeatFluxIC
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFluxIC::~wallHeatFluxIC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFluxIC::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    resetName(typeName);
    resetLocalObjectName(typeName);

    return true;
}


bool Foam::functionObjects::wallHeatFluxIC::execute()
{

    tmp<volScalarField> tT;
    if (mesh_.foundObject<volScalarField>("T"))
    {
        tT = mesh_.lookupObject<volScalarField>("T");
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find temperature field in the "
            << "database" << exit(FatalError);
    }
    const volScalarField& T = tT();
    word name(type());
    tmp<volVectorField> tgradT = fvc::grad(T);
    return store(name, calcHeatFlux(tgradT()));
}


bool Foam::functionObjects::wallHeatFluxIC::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& wallHeatFluxIC =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFluxIC.boundaryField()[patchi];

        scalar minHfp = gMin(hfp);
        scalar maxHfp = gMax(hfp);

        if (Pstream::master())
        {
            file() << mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
