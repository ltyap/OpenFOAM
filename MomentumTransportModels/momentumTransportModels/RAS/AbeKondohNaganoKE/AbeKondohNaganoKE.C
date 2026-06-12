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

#include "AbeKondohNaganoKE.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE<BasicMomentumTransportModel>::yStar() const
{
    return pow(this->nu()*epsilon_, 0.25)*y_/this->nu();
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE<BasicMomentumTransportModel>::Rt() const
{
    return sqr(k_)/(this->nu()*epsilon_);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE<BasicMomentumTransportModel>::fMu() const
{
    tmp<volScalarField> fMuPlus = sqr(scalar(1) - exp(-yStar()/14.0))
                                *(
                                    scalar(1) + (5.0/(pow(Rt(), 0.75) + small))
                                    *exp(-sqr(Rt()/200.0))
                                );
    return min(scalar(1.0), fMuPlus);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE<BasicMomentumTransportModel>::fEps2() const
{

    tmp<volScalarField> RtPlus = min(sqr(Rt()/6.5), scalar(30.0));
    return (scalar(1)-0.3*exp(-RtPlus))*sqr(scalar(1)-exp(-yStar()/3.1));
}

template<class BasicMomentumTransportModel>
void AbeKondohNaganoKE<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    // Info << "k" << endl;
    // Info << this->k_ << endl;    
    // Info << "epsilon" << endl;
    // Info << this->epsilon_ << endl;
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
AbeKondohNaganoKE<BasicMomentumTransportModel>::AbeKondohNaganoKE
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
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            1.5
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9
        )
    ),
    Ceps3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.4
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.4
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    y_(wallDist::New(this->mesh_).y())
    
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool AbeKondohNaganoKE<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Ceps3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void AbeKondohNaganoKE<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }
 
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
    // const volScalarField fEps2(this->fEps2());

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );    
    tgradU.clear();
    // Info << "Before solving kTheta and epsTheta: " << endl;
    // Info << "G" << endl;
    // Info << G << endl;
    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Info << "After BC update: " << endl;
    // Info << "G" << endl;
    // Info << G << endl;

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*Ceps1_ - Ceps3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(Ceps2_*fEps2()*alpha*rho*epsilon_/k_, epsilon_)
      + fvModels.source(alpha, rho, epsilon_)
    );
    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);
    
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
