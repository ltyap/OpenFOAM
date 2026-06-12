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

#include "AKN_DWX.H"
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
tmp<volScalarField> AKN_DWX<BasicMomentumTransportModel>::fd1() const
{
    return sqr(1 - exp(-this->yStar()/1.7));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_DWX<BasicMomentumTransportModel>::fd2() const
{
    tmp<volScalarField> f1{scalar(1) - 0.3*exp(-sqr(this->Rt()/6.5))};
    tmp<volScalarField> f2{scalar(1) - exp(-this->yStar()/5.8)};
    return (1/Cd2_)*(this->Ceps2_*f1 - scalar(1))*sqr(f2);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_DWX<BasicMomentumTransportModel>::R() const
{
    tmp<volScalarField> tauU = this->k()/this->epsilon();
    tmp<volScalarField> tauT = kTheta_/epsilonTheta_;

    return tauT/(tauU + dimensionedScalar(dimTime, rootVSmall));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AKN_DWX<BasicMomentumTransportModel>::fAlpha
() const
{
    tmp<volScalarField> fLambda = sqr(scalar(1) - exp(-this->yStar()/16.0))
                                    *(1 + 3/(pow(this->Rt(), 0.75) + small));
    return fLambda*2.0*sqrt(R());
}

template<class BasicMomentumTransportModel>
void AKN_DWX<BasicMomentumTransportModel>::correctAlphat()
{
    Dt_ = Clambda_*fAlpha()*sqr(this->k())/this->epsilon();
    Dt_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(Dt_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
AKN_DWX<BasicMomentumTransportModel>::AKN_DWX
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
    AbeKondohNaganoKE<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
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
            1.654
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
    Cd2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd2",
            this->coeffDict_,
            0.9
        )
    ),
    sigmakTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmakTheta",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEpsTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEpsTheta",
            this->coeffDict_,
            1.0
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ) ,

    D_ (this->mesh_.objectRegistry::template lookupObject<volScalarField>("alpha"))
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
bool AKN_DWX<BasicMomentumTransportModel>::read()
{
    if (AbeKondohNaganoKE<BasicMomentumTransportModel>::read())
    {
        Cp1_.readIfPresent(this->coeffDict());
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());
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
void AKN_DWX<BasicMomentumTransportModel>::correctThermal()
{

    const surfaceScalarField& phi = this->phi_;
    const volScalarField& T = this->mesh_.objectRegistry::template lookupObject<volScalarField>("T");
    const volScalarField& k = this->k_; 
    const volScalarField& epsilon = this->epsilon_; 

    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    
    tmp<volVectorField> tgradT = fvc::grad(T);
    volScalarField::Internal GTheta
    (
        "GTheta",
        Dt_.v()*(tgradT().v() & tgradT().v())
    );    
    tgradT.clear();

    // Update epsilon and GTheta at the wall
    // Info << GTheta << endl;
    epsilonTheta_.boundaryFieldRef().updateCoeffs();
    // Info << GTheta << endl;
    // Info << "after update of epsilonTheta" << endl;
    
    tmp<fvScalarMatrix> epsThetaEqn
    (
        fvm::ddt(epsilonTheta_)
      + fvm::div(phi, epsilonTheta_)
      - fvm::laplacian(DepsilonThetaEff(), epsilonTheta_)
     ==
        Cp1_*GTheta*sqrt((epsilonTheta_*epsilon())/(kTheta_*k()))
      - fvm::Sp(Cd1_*fd1()*epsilonTheta_/kTheta_, epsilonTheta_)
      - fvm::Sp(Cd2_*fd2()*epsilon/k, epsilonTheta_)
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

    AKN_DWX<BasicMomentumTransportModel>::correctAlphat();
}

template<class BasicMomentumTransportModel>
void AKN_DWX<BasicMomentumTransportModel>::correct()
{
    AbeKondohNaganoKE<BasicMomentumTransportModel>::correct();
    correctThermal();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
