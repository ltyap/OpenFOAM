/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation, 2020 Timofey Mukha
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

#include "channelIndex.H"
#include "boolList.H"
#include "syncTools.H"
#include "OFstream.H"
#include "meshTools.H"
#include "Time.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::vector::components,
        3
    >::names[] =
    {
        "x",
        "y",
        "z"
    };
}
 
const Foam::NamedEnum<Foam::vector::components, 3>
     Foam::channelIndex::vectorComponentsNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Determines face blocking
void Foam::channelIndex::walkOppositeFaces
(
    const polyMesh& mesh,
    const labelList& startFaces,
    boolList& blockedFace
)
{
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const label nBnd = mesh.nFaces() - mesh.nInternalFaces();

    DynamicList<label> frontFaces(startFaces);
    forAll(frontFaces, i)
    {
        label facei = frontFaces[i];
        blockedFace[facei] = true;
    }

    while (returnReduce(frontFaces.size(), sumOp<label>()) > 0)
    {
        // Transfer across.
        boolList isFrontBndFace(nBnd, false);
        forAll(frontFaces, i)
        {
            label facei = frontFaces[i];

            if (!mesh.isInternalFace(facei))
            {
                isFrontBndFace[facei-mesh.nInternalFaces()] = true;
            }
        }
        syncTools::swapBoundaryFaceList(mesh, isFrontBndFace);

        // Add
        forAll(isFrontBndFace, i)
        {
            label facei = mesh.nInternalFaces()+i;
            if (isFrontBndFace[i] && !blockedFace[facei])
            {
                blockedFace[facei] = true;
                frontFaces.append(facei);
            }
        }

        // Transfer across cells
        DynamicList<label> newFrontFaces(frontFaces.size());

        forAll(frontFaces, i)
        {
            label facei = frontFaces[i];

            {
                const cell& ownCell = cells[mesh.faceOwner()[facei]];

                label oppositeFacei = ownCell.opposingFaceLabel(facei, faces);

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Face:" << facei << " owner cell:" << ownCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFacei])
                    {
                        blockedFace[oppositeFacei] = true;
                        newFrontFaces.append(oppositeFacei);
                    }
                }
            }

            if (mesh.isInternalFace(facei))
            {
                const cell& neiCell = mesh.cells()[mesh.faceNeighbour()[facei]];

                label oppositeFacei = neiCell.opposingFaceLabel(facei, faces);

                if (oppositeFacei == -1)
                {
                    FatalErrorInFunction
                        << "Face:" << facei << " neighbour cell:" << neiCell
                        << " is not a hex?" << abort(FatalError);
                }
                else
                {
                    if (!blockedFace[oppositeFacei])
                    {
                        blockedFace[oppositeFacei] = true;
                        newFrontFaces.append(oppositeFacei);
                    }
                }
            }
        }

        frontFaces.transfer(newFrontFaces);
    }
}


// Calculate regions.
void Foam::channelIndex::calcLayeredRegions
(
    const polyMesh& mesh,
    const boolList& blockedFace
)
{


    if (false)
    {
        OFstream str(mesh.time().path()/"blockedFaces.obj");
        label vertI = 0;
        forAll(blockedFace, facei)
        {
            if (blockedFace[facei])
            {
                const face& f = mesh.faces()[facei];
                forAll(f, fp)
                {
                    meshTools::writeOBJ(str, mesh.points()[f[fp]]);
                }
                str<< 'f';
                forAll(f, fp)
                {
                    str << ' ' << vertI+fp+1;
                }
                str << nl;
                vertI += f.size();
            }
        }
    }


    // Do analysis for connected regions
    cellRegion_.reset(new regionSplit(mesh, blockedFace));

    Info<< "Detected " << cellRegion_().nRegions() << " layers." << nl << endl;

    // Sum number of entries per region
    regionCount_ = regionSum(scalarField(mesh.nCells(), 1.0));

    tmp<scalarField> radiustmp{calcCellRadius(mesh.cellCentres())};
    scalarField& radius = radiustmp.ref();

    // Average cell radius to determine ordering.
    scalarField regionRadius
    {
        regionSum(radius)
      / regionCount_
    };

    SortableList<scalar> sortComponent{regionRadius};

    sortMap_ = sortComponent.indices();
    radiusInternal_ = sortComponent;
}


void Foam::channelIndex::findBottomPatchIndices
(
    const polyMesh& mesh,
    const wordList& patchNames
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    for (word i : patchNames)
    {
        const label patchI = bMesh.findPatchID(i);

        if (patchI == -1)
        {
            FatalErrorInFunction
                << "Illegal patch " << i
                << ". Valid patches are " << bMesh.name()
                << exit(FatalError);
        }

        bottomPatchIndices_.append(patchI);
    }
}

void Foam::channelIndex::findBottomPatchIndices
(
    const polyMesh& mesh,
    const labelList& startFaces
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    for (label i : startFaces)
    {
        const label patchI = bMesh.whichPatch(i);

        if (Foam::findIndex(bottomPatchIndices_, patchI ) == -1)
        {
            bottomPatchIndices_.append(patchI);
        }

    }
}


void Foam::channelIndex::findTopPatchIndices
(
    const polyMesh& mesh,
    const boolList& blockedFace
)
{
    const polyBoundaryMesh & bMesh = mesh.boundaryMesh();
    const label nInt =  mesh.nInternalFaces();
    const label nBnd = mesh.nFaces() - nInt;

    for (label i=0; i<nBnd; i++)
    {
        label faceI = nInt + i;

        if (blockedFace[faceI])
        {
            label patchI = bMesh.whichPatch(faceI);
           
            if (
                    (Foam::findIndex(bottomPatchIndices_, patchI) == -1)
                    &&
                    (Foam::findIndex(topPatchIndices_, patchI) == -1)
                )
            {
                topPatchIndices_.append(patchI); 
            }

        }
    }

    if (topPatchIndices_.size() == 0)
    {
        FatalErrorInFunction
            << "Could not find the top patch(es)."
            << exit(FatalError);
    }

    Info<< "The top patches are: ";
    wordList patchNames = bMesh.names();
    for (label i : topPatchIndices_)
    {
        Info<< patchNames[i] << " ";
    }
    Info << nl;
}

void Foam::channelIndex::checkPatchSizes
(
    const polyBoundaryMesh& bMesh
)
{
    label nBottom = 0;
    label nTop = 0;

    for (label i : bottomPatchIndices_)
    {
        nBottom += bMesh[i].size();
    }
    for (label i : topPatchIndices_)
    {
        nTop += bMesh[i].size();
    }
    
    if (nBottom != nTop)
    {
        FatalErrorInFunction
            << "The total number of faces on the bottom and top patches are"
            << "unequal. Bottom faces: "<< nBottom << " Top faces: "
            << nTop << "."
            << exit(FatalError);
    }
}


Foam::tmp<Foam::scalarField>
Foam::channelIndex::y(const polyBoundaryMesh& bMesh) const
{
    tmp<scalarField> ytmp(new scalarField(radiusInternal_.size() + 2, Foam::Zero));
    scalarField& y = ytmp.ref();

    for (label i=0; i<radiusInternal_.size(); ++i)
    {
        y[i+1] = radiusInternal_[i];
    }

    // We rely on all faces of the corresponding wall boundary to 
    // have the same radius. Reasonable for annular channel flow.
    vectorField fc = bMesh[bottomPatchIndices_[0]].faceCentres();

    // y[0] = Foam::sqrt(Foam::sqr(fc[0].component(dir1_)) + Foam::sqr(fc[0].component(dir2_)));
    y[0] = radius(fc[0]);
    fc = bMesh[topPatchIndices_[0]].faceCentres();
    // y[y.size()-1] = Foam::sqrt(Foam::sqr(fc[0].component(dir1_)) + Foam::sqr(fc[0].component(dir2_)));
    y[y.size()-1] = radius(fc[0]);

    const scalar innerWallRadius{y[0]};
    y -= innerWallRadius;

    return ytmp;   
}

Foam::tmp<Foam::scalarField> 
Foam::channelIndex::calcCellRadius(const vectorField& cellField) const
{
    tmp<scalarField> cellRadiustmp{new scalarField (cellField.size(), Foam::Zero)};
    scalarField& cellRadius = cellRadiustmp.ref();
    // cellRadius = Foam::sqrt(Foam::sqr(cellField.component(dir1_)) + Foam::sqr(cellField.component(dir2_)));
    cellRadius = this->radius(cellField);

    return cellRadiustmp;
}

Foam::scalarField Foam::channelIndex::radius(const vectorField& cellField) const
{
    return Foam::sqrt(Foam::sqr(cellField.component(dir1_)) + Foam::sqr(cellField.component(dir2_)));
}

Foam::scalar Foam::channelIndex::radius(const vector& point) const
{
    return Foam::sqrt(Foam::sqr(point.component(dir1_)) + Foam::sqr(point.component(dir2_)));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    axialDir_(vectorComponentsNames_.read(dict.lookup("axis"))),
    dir1_(1),
    dir2_(2)
{

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const wordList patchNames(dict.lookup("patches"));

    // Get the seed patch indices from the patch names
    findBottomPatchIndices(mesh, patchNames);

    // Sum the number of faces on the seed patches
    label nFaces = 0;

    forAll(patchNames, i)
    {
        nFaces += bMesh[bottomPatchIndices_[i]].size();
    }

    labelList startFaces(nFaces);
    nFaces = 0;

    forAll(patchNames, i)
    {
        const polyPatch& pp = bMesh[patchNames[i]];

        forAll(pp, j)
        {
            startFaces[nFaces++] = pp.start()+j;
        }
    }

    boolList blockedFace(mesh.nFaces(), false);
    walkOppositeFaces
    (
        mesh,
        startFaces,
        blockedFace
    );

    // Find 
    findTopPatchIndices(mesh, blockedFace);

    checkPatchSizes(bMesh);

    // Calculate regions.
    calcLayeredRegions(mesh, blockedFace);
}

/*
Foam::channelIndex::channelIndex
(
    const polyMesh& mesh,
    const labelList& startFaces,
    //const bool symmetric,
    const direction dir
)
:
    //symmetric_(symmetric),
    dir_(dir)
{
    boolList blockedFace(mesh.nFaces(), false);
    walkOppositeFaces
    (
        mesh,
        startFaces,
        blockedFace
    );
    

    findBottomPatchIndices(mesh, startFaces);
    findTopPatchIndices(mesh, blockedFace);

    checkPatchSizes(mesh.boundaryMesh());
    
    // Calculate regions.
    calcLayeredRegions(mesh, blockedFace);
}
*/

// ************************************************************************* //
