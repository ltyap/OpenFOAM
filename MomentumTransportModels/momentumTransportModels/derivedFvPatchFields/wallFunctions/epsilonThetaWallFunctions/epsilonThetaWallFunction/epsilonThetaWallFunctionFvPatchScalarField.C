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

\*---------------------------------------------------------------------------*/

#include "epsilonThetaWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonThetaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::epsilonThetaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilonTheta =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilonTheta.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonThetaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonThetaWallFunctionFvPatchScalarField& epf = epsilonThetaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilonTheta =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilonTheta.boundaryField();

    const fvMesh& mesh = epsilonTheta.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonThetaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchi = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    // G_.setSize(internalField().size(), 0.0);
    GTheta_.setSize(internalField().size(), 0.0);
    epsilonTheta_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::epsilonThetaWallFunctionFvPatchScalarField&
Foam::epsilonThetaWallFunctionFvPatchScalarField::epsilonThetaPatch(const label patchi)
{
    const volScalarField& epsilonTheta =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilonTheta.boundaryField();

    const epsilonThetaWallFunctionFvPatchScalarField& epf =
        refCast<const epsilonThetaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonThetaWallFunctionFvPatchScalarField&>(epf);
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const momentumTransportModel& turbulence,
    // scalarField& G0,
    scalarField& GTheta0,
    scalarField& epsilonTheta0
)
{
    // Accumulate all of the G and epsilonTheta contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonThetaWallFunctionFvPatchScalarField& epf = epsilonThetaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            // epf.calculate(turbulence, w, epf.patch(), G0, GTheta0, epsilonTheta0);
            epf.calculate(turbulence, w, epf.patch(), GTheta0, epsilonTheta0);
        }
    }

    // Apply zero-gradient condition for epsilonTheta
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonThetaWallFunctionFvPatchScalarField& epf = epsilonThetaPatch(patchi);

            epf == scalarField(epsilonTheta0, epf.patch().faceCells());
        }
    }
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::calculate
(
    const momentumTransportModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    // scalarField& G0,
    scalarField& GTheta0,
    scalarField& epsilonTheta0
)
{
    const label patchi = patch.index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    // const tmp<volScalarField> tk = turbModel.k();
    // const volScalarField& k = tk();

    const volScalarField& kTheta = db().lookupObject<volScalarField>("kTheta");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const volScalarField& alpha = db().lookupObject<volScalarField>("alpha");
    const scalarField& alphaw = alpha.boundaryField()[patchi];

    const volScalarField& alphat = db().lookupObject<volScalarField>("alphat");
    const scalarField& alphatw = alphat.boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu75 = pow(nutw.Cmu(), 0.75);

    // Set epsilonTheta and G
    forAll(alphatw, facei)
    {
        const label celli = patch.faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(kTheta[celli])/nuw[facei];

        const scalar w = cornerWeights[facei];

        // The first part of the algorithm is most likely wrong but we dont need this for now......
        if (yPlus > nutw.yPlusLam())
        {
            epsilonTheta0[celli] +=
                w*Cmu75*pow(kTheta[celli], 1.5)/(nutw.kappa()*y[facei]);
            
            // G0[celli] +=
            //     w
            //    *(nutw[facei] + nuw[facei])
            //    *magGradUw[facei]
            //    *Cmu25*sqrt(k[celli])
            //   /(nutw.kappa()*y[facei]);
            
            GTheta0[celli] +=
                w
               *(alphatw[facei] + alphaw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(kTheta[celli])
              /(nutw.kappa()*y[facei]);
        }
        else
        {
            // Info << alphaw << endl;
            // Info << kTheta[celli] << endl;
            // Info << "-------------" << endl;
            // Info << alphaw[facei] << endl;
            // Info << "-------------" << endl;
            // Info << y[facei] << endl;
            // Info << "-------------" << endl;

            epsilonTheta0[celli] += w*2.0*kTheta[celli]*alphaw[facei]/sqr(y[facei]);
            // Info << epsilonTheta0[celli] << endl;
            // Info << "***********************" << endl;

        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonThetaWallFunctionFvPatchScalarField::
epsilonThetaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    // G_(),
    GTheta_(),
    epsilonTheta_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


Foam::epsilonThetaWallFunctionFvPatchScalarField::
epsilonThetaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    // G_(),
    GTheta_(),
    epsilonTheta_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::epsilonThetaWallFunctionFvPatchScalarField::
epsilonThetaWallFunctionFvPatchScalarField
(
    const epsilonThetaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    // G_(),
    GTheta_(),
    epsilonTheta_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


Foam::epsilonThetaWallFunctionFvPatchScalarField::
epsilonThetaWallFunctionFvPatchScalarField
(
    const epsilonThetaWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    // G_(),
    GTheta_(),
    epsilonTheta_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Foam::scalarField& Foam::epsilonThetaWallFunctionFvPatchScalarField::G(bool init)
// {
//     if (patch().index() == master_)
//     {
//         if (init)
//         {
//             G_ = 0.0;
//         }

//         return G_;
//     }

//     return epsilonThetaPatch(master_).G();
// }

Foam::scalarField& Foam::epsilonThetaWallFunctionFvPatchScalarField::GTheta(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            GTheta_ = 0.0;
        }

        return GTheta_;
    }

    return epsilonThetaPatch(master_).GTheta();
}


Foam::scalarField& Foam::epsilonThetaWallFunctionFvPatchScalarField::epsilonTheta
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilonTheta_ = 0.0;
        }

        return epsilonTheta_;
    }

    return epsilonThetaPatch(master_).epsilonTheta(init);
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        // calculateTurbulenceFields(turbModel, G(true), GTheta(true), epsilonTheta(true));
        calculateTurbulenceFields(turbModel, GTheta(true), epsilonTheta(true));

    }

    // const scalarField& G0 = this->G();
    const scalarField& GTheta0 = this->GTheta();
    const scalarField& epsilonTheta0 = this->epsilonTheta();

    typedef DimensionedField<scalar, volMesh> FieldType;
    
    // FieldType& G =
    //     const_cast<FieldType&>
    //     (
    //         db().lookupObject<FieldType>(turbModel.GName())
    //     );

    FieldType& GTheta =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>("GTheta")
        );

    FieldType& epsilonTheta = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];
        
        // G[celli] = G0[celli];
        GTheta[celli] = GTheta0[celli];
        epsilonTheta[celli] = epsilonTheta0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        // calculateTurbulenceFields(turbModel, G(true), GTheta(true), epsilonTheta(true));
        calculateTurbulenceFields(turbModel, GTheta(true), epsilonTheta(true));
    }

    // const scalarField& G0 = this->G();
    const scalarField& GTheta0 = this->GTheta();
    const scalarField& epsilonTheta0 = this->epsilonTheta();

    typedef DimensionedField<scalar, volMesh> FieldType;

    // FieldType& G =
    //     const_cast<FieldType&>
    //     (
    //         db().lookupObject<FieldType>(turbModel.GName())
    //     );

    FieldType& GTheta =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>("GTheta")
        );

    FieldType& epsilonTheta = const_cast<FieldType&>(internalField());

    scalarField& epsilonThetaf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            // G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            GTheta[celli] = (1.0 - w)*GTheta[celli] + w*GTheta0[celli];
            epsilonTheta[celli] = (1.0 - w)*epsilonTheta[celli] + w*epsilonTheta0[celli];
            epsilonThetaf[facei] = epsilonTheta[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonThetaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintEpsilon(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilonTheta
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintEpsilon.append(epsilonTheta[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintEpsilon)
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonThetaWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
