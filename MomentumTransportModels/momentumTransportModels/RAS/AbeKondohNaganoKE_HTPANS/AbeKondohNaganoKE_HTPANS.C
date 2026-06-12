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

#include "AbeKondohNaganoKE_HTPANS.H"
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
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::yStar() const
{
    return pow(this->nu()*epsilon_, 0.25)*y_/this->nu();
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::Rt() const
{
    return sqr(k_)/(this->nu()*epsilon_);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::fMu() const
{
    tmp<volScalarField> fMuPlus = sqr(scalar(1) - exp(-yStar()/14.0))
                                *(
                                    scalar(1) + (5.0/(pow(Rt(), 0.75) + small))
                                    *exp(-sqr(Rt()/200.0))
                                );
    return min(scalar(1.0), fMuPlus);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::fEps2() const
{

    tmp<volScalarField> RtPlus = min(sqr(Rt()/6.5), scalar(30.0));
    return (scalar(1)-0.3*exp(-RtPlus))*sqr(scalar(1)-exp(-yStar()/3.1));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::fp1() const
{
    return sqr(1 - exp(-this->yStar()));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::fd1() const
{
    return sqr(1 - exp(-this->yStar()));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::fd2() const
{
    tmp<volScalarField> f1{scalar(1) - 0.3*exp(-sqr(Rt()/6.5))};
    tmp<volScalarField> f2{scalar(1) - exp(-yStar()/5.7)};
    return (1/Cd2_)*(Ceps2_*f1 - scalar(1))*sqr(f2);

}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::R() const
{
    tmp<volScalarField> tauU = k_/epsilon_;
    tmp<volScalarField> tauT = kTheta_/epsilonTheta_;

    return tauT/(tauU + dimensionedScalar(dimTime, rootVSmall));
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::fAlpha
() const
{

    tmp<volScalarField> fLambda = (scalar(1) - exp(-yStar()/14.0))*(scalar(1) - exp(-yStar()*sqrt(Pr_)/14.0));
    tmp<volScalarField> fTheta2 = exp(-sqr(Rt()/200.0));

    tmp<volScalarField> cTheta1 = 2*R()/(Cgamma_ + R());
    tmp<volScalarField> cTheta2 = 3.0*sqrt(2*R())/(Pr_*pow(Rt(), 0.75) + small);

    return fLambda*(cTheta1 + cTheta2*fTheta2);
}

template<class BasicMomentumTransportModel>
void AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);

    Dt_ = Clambda_*fAlpha()*sqr(k_)/epsilon_;
    Dt_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(Dt_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::AbeKondohNaganoKE_HTPANS
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
            0.95
        )
    ),
    Cp2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp2",
            this->coeffDict_,
            0.6
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
    Cgamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cgamma",
            this->coeffDict_,
            0.5
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

    Pr_    
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Pr",
            transportProperties_,
            0.71
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

    y_(wallDist::New(this->mesh_).y()),

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

    D_ (this->mesh_.objectRegistry::template lookupObject<volScalarField>("alpha")),

    fK_
    (
        IOobject
        (
            IOobject::groupName("fK", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ),

    fEps_
    (
        IOobject
        (
            IOobject::groupName("fEps", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, 1.0)
    ), 

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
        dimensionedScalar(dimless, 1.0)
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
        dimensionedScalar(dimless, 1.0)
    )  

{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(kTheta_, kThetaMin_);
    bound(epsilonTheta_, epsilonThetaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        Clambda_.readIfPresent(this->coeffDict());
        Cp1_.readIfPresent(this->coeffDict());
        Cp2_.readIfPresent(this->coeffDict());
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());
        Cgamma_.readIfPresent(this->coeffDict());
        sigmakTheta_.readIfPresent(this->coeffDict());
        sigmaEpsTheta_.readIfPresent(this->coeffDict());
        Pr_.readIfPresent(transportProperties_);

        return true;
    }
    else
    {
        return false;
    }
}
template<class BasicMomentumTransportModel>
void AbeKondohNaganoKE_HTPANS<BasicMomentumTransportModel>::correct()
{

    if (!this->turbulence_)
    {
        return;
    }
 
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const surfaceScalarField& phi = this->phi_;
    const volVectorField& U = this->U_;
    const volScalarField& T = this->mesh_.objectRegistry::template lookupObject<volScalarField>("T");
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();
    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );    
    tgradU.clear();
 
    epsilon_.boundaryFieldRef().updateCoeffs();

    const volScalarField Ceps2Star_ =  Ceps1_ + (fK_/fEps_)*(Ceps2_*fEps2() - Ceps1_);

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*Ceps1_ - Ceps3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(Ceps2Star_*alpha()*rho()*epsilon_()/k_(), epsilon_)
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


    volScalarField GammaEff("GammaEff", D_ + Dt_);

    const volScalarField fp1(this->fp1());
    const volScalarField fd1(this->fd1());
    const volScalarField fd2(this->fd2());


    tmp<volVectorField> tgradT = fvc::grad(T);

    volScalarField::Internal GTheta
    (
        "GTheta",
        Dt_.v()*(tgradT().v() & tgradT().v())
    );    
    tgradT.clear();

    // Update epsilon and GTheta at the wall
    epsilonTheta_.boundaryFieldRef().updateCoeffs();

    // Info << "After: "<< endl;
    // Info << "G:" << endl;
    // Info << G << endl;
    // Info << "GTheta:" << endl;
    // Info << GTheta << endl;
    // Info << "epsilon" << endl;
    // Info << epsilon << endl;
    // Info << "***************" << endl;

    const volScalarField Cd1Star_ =  Cp1_*fp1 + (fKTheta_/fEpsTheta_)*(Cd1_*fd1 - Cp1_*fp1);
    const volScalarField Cd2Star_ =  Cp2_ + (fK_/fEps_)*(Cd2_*fd2 - Cp2_);


    tmp<fvScalarMatrix> epsThetaEqn
    (
        fvm::ddt(epsilonTheta_)
      + fvm::div(phi, epsilonTheta_)
      - fvm::laplacian(DepsilonThetaEff(), epsilonTheta_)
     ==
        Cp1_*fp1()*GTheta*epsilonTheta_/kTheta_
      + Cp2_*G*epsilonTheta_/k_
      - fvm::Sp(Cd1Star_*epsilonTheta_()/kTheta_(), epsilonTheta_)
      - fvm::Sp(Cd2Star_*epsilon_()/k_(), epsilonTheta_)
    );

    Info << Pr_<< endl;

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

    correctNut();

}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
