/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "AKN_MMPANS.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_MMPANS<BasicMomentumTransportModel>::fd2() const
{
    tmp<volScalarField> f1{scalar(1) - 0.3*exp(-sqr(this->Rt()/6.5))};
    tmp<volScalarField> f2{scalar(1) - exp(-this->yStar()/5.7)};
    return (this->Ceps2_*f1 - scalar(1))*sqr(f2);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_MMPANS<BasicMomentumTransportModel>::f1Theta() const
{
    return (scalar(1) - exp(-this->yStar()/14.0))
          *(scalar(1) - exp(-this->yStar()*sqrt(Pr_)/19.0));
}
template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_MMPANS<BasicMomentumTransportModel>::R() const
{
    tmp<volScalarField> tauU = this->k()/this->epsilon();
    tmp<volScalarField> tauT = kTheta_/epsilonTheta_;
    return tauT/(tauU + dimensionedScalar(dimTime, rootVSmall));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_MMPANS<BasicMomentumTransportModel>::fAlpha() const
{
    return f1Theta()
       *(
            Binf_ 
            + exp(-sqr(this->Rt()/500))*(2*R()/(R() + Cgamma_))
            + exp(-sqr(this->Rt()/200))*sqrt(2*R())*(1.3/(Pr_*pow(this->Rt(), 0.75) + small))
        );
}

template<class BasicMomentumTransportModel>
void AKN_MMPANS<BasicMomentumTransportModel>::correctAlphat()
{

    Dt_ = Clambda_*fAlpha()*sqr(this->k())/this->epsilon();
    Dt_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(Dt_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
AKN_MMPANS<BasicMomentumTransportModel>::AKN_MMPANS
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    AKNPANS<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",    // dictionary name
            this->runTime_.constant(),     // dict is found in "constant"
            this->mesh_,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            this->coeffDict_,
            0.10
        )
    ),
    Cp1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp1",
            this->coeffDict_,
            0.925
        )
    ),
    Cp2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp2",
            this->coeffDict_,
            0.9
        )
    ),
    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            this->coeffDict_,
            1.0
        )
    ),
    Cgamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cgamma",
            this->coeffDict_,
            0.3
        )
    ),
    Binf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Binf",
            this->coeffDict_,
            0.9
        )
    ),
    Pr_    
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Pr",
            transportProperties_,
            0.71
        )
    ),
    sigmakTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmakTheta",
            this->coeffDict_,
            1.6
        )
    ),
    sigmaEpsTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsTheta",
            this->coeffDict_,
            1.6
        )
    ),
    kThetaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kThetaMin",
            this->RASDict_,
            sqr(dimTemperature),
            small
        )
    ),

    epsilonThetaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilonThetaMin",
            this->RASDict_,
            kThetaMin_.dimensions()/dimTime,
            small
        )
    ),
    kTheta_
    (
        IOobject
        (
            IOobject::groupName("kTheta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilonTheta_
    (
        IOobject
        (
            IOobject::groupName("epsilonTheta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Dt_
    (
        IOobject
        (
            IOobject::groupName("alphat", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar ("alphat", dimViscosity, SMALL)
    ),

    D_ (this->mesh_.objectRegistry::template lookupObject<volScalarField>("alpha")),

    fKTheta_
    (
        IOobject
        (
            IOobject::groupName("fKTheta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("fKTheta", dimless, 1.0)
    ),

    fEpsTheta_
    (
        IOobject
        (
            IOobject::groupName("fEpsTheta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("fEpsTheta", dimless, 1.0)
    )  
{
    bound(kTheta_, kThetaMin_);
    bound(epsilonTheta_, epsilonThetaMin_);
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool AKN_MMPANS<BasicMomentumTransportModel>::read()
{
    if (AKNPANS<BasicMomentumTransportModel>::read())
    {
        Clambda_.readIfPresent(this->coeffDict());
        Cp1_.readIfPresent(this->coeffDict());
        Cp2_.readIfPresent(this->coeffDict());
        Cd1_.readIfPresent(this->coeffDict());
        Cgamma_.readIfPresent(this->coeffDict());
        Binf_.readIfPresent(this->coeffDict());
        sigmakTheta_.readIfPresent(this->coeffDict());
        sigmaEpsTheta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}
template<class BasicMomentumTransportModel>
void AKN_MMPANS<BasicMomentumTransportModel>::correctThermal()
{

    // Local references
    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    const volScalarField& T = this->mesh_.objectRegistry::template lookupObject<volScalarField>("T");
    volScalarField& nut = this->nut_;
    const volScalarField k(this->k()); 
    const volScalarField epsilon(this->epsilon()); 

    // const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );
   
    tmp<volTensorField> tgradU = fvc::grad(U);
    
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );    
    tgradU.clear();

    tmp<volVectorField> tgradT = fvc::grad(T);
    volScalarField::Internal GTheta
    (
        "GTheta",
        Dt_.v()*(tgradT().v() & tgradT().v())
    );    
    tgradT.clear();

    // Update epsilon and G at the wall
    epsilonTheta_.boundaryFieldRef().updateCoeffs();


    const volScalarField Cd1Star_ =  Cp1_ + (fKTheta_/fEpsTheta_)*(Cd1_ - Cp1_);
    const volScalarField Cd2Star_ =  Cp2_ + (this->fK_/this->fEps_)*(fd2() - Cp2_);

    tmp<fvScalarMatrix> epsThetaEqn
    (
        fvm::ddt(epsilonTheta_)
      + fvm::div(phi, epsilonTheta_)
      - fvm::laplacian(DepsilonThetaEff(), epsilonTheta_)
     ==
        Cp1_*GTheta*epsilonTheta_()/kTheta_()
      + Cp2_*G*epsilonTheta_()/k()
      - fvm::Sp(Cd1Star_()*epsilonTheta_()/kTheta_(), epsilonTheta_)
      - fvm::Sp(Cd2Star_()*epsilon()/k(), epsilonTheta_)
    );
    epsThetaEqn.ref().relax();
    fvConstraints.constrain(epsThetaEqn.ref());
    epsThetaEqn.ref().boundaryManipulate(epsilonTheta_.boundaryFieldRef());
    solve(epsThetaEqn);
    fvConstraints.constrain(epsilonTheta_);
    bound(epsilonTheta_, epsilonThetaMin_);
    
    tmp<fvScalarMatrix> kThetaEqn
    (
        fvm::ddt(kTheta_)
      + fvm::div(phi, kTheta_)
      - fvm::laplacian(DkThetaEff(), kTheta_)
     ==
        GTheta
      - fvm::Sp(epsilonTheta_()/kTheta_(), kTheta_)
    );

    kThetaEqn.ref().relax();
    fvConstraints.constrain(kThetaEqn.ref());
    solve(kThetaEqn);
    fvConstraints.constrain(kTheta_);
    bound(kTheta_, kThetaMin_);

    correctAlphat();
}

template<class BasicMomentumTransportModel>
void AKN_MMPANS<BasicMomentumTransportModel>::correct()
{
    AKNPANS<BasicMomentumTransportModel>::correct();
    correctThermal();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
