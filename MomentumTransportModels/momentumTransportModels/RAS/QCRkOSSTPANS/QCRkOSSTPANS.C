/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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

#include "QCRkOSSTPANS.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
void QCRkOSSTPANS<BasicMomentumTransportModel>::correctNonlinearStress
(const volTensorField& gradU)
{
    volSymmTensorField R(-eddyViscosity<RASModel<BasicMomentumTransportModel>>::sigma());
    volTensorField W(skew(gradU));
    volScalarField magGradU(mag(gradU));
    volTensorField Wbar(2.0*W/(magGradU + dimensionedScalar("small", magGradU.dimensions(), ROOTVSMALL)));

    this->nonlinearStress_ = Cr1_*twoSymm(Wbar&R);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
QCRkOSSTPANS<BasicMomentumTransportModel>::QCRkOSSTPANS
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
    Foam::kOmegaSSTPANS
    <
        nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),

    Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            this->coeffDict_,
            0.3
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicMomentumTransportModel>
bool QCRkOSSTPANS<BasicMomentumTransportModel>::read()
{
    if (Foam::kOmegaSSTPANS
    <
        nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >::read())
    {
        Cr1_.readIfPresent(this->coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void QCRkOSSTPANS<BasicMomentumTransportModel>::correct()
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

    BasicMomentumTransportModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();

    volScalarField S2(2*magSqr(symm(tgradU())));

    volScalarField::Internal G(
        this->GName(),
        (nut()*dev(twoSymm(tgradU()())) - this->nonlinearStress_) && tgradU()()
    );
    

    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();

    volScalarField::Internal GbyNu(G/nut());


    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_*(this->fOmega_/this->fK_))
        *(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F2(this->F2());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));
        volScalarField::Internal betaPrime((gamma*this->betaStar_) 
                                        - (gamma*this->betaStar_/this->fOmega_) 
                                         + (beta/this->fOmega_)
                                            );        
        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, this->omega_)
          + fvm::div(alphaRhoPhi, this->omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNu,
                (this->c1_/this->a1_)*this->betaStar_*this->omega_()
               *max(this->a1_*this->omega_(), this->b1_*F2()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
          - fvm::Sp(alpha()*rho()*betaPrime*this->omega_(), this->omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvModels.source(alpha, rho, this->omega_)
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, this->k_)
      + fvm::div(alphaRhoPhi, this->k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_)
     ==
        alpha()*rho()*this->Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
      - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, F2), this->k_)
      + this->kSource()
      + fvModels.source(alpha, rho, this->k_)
    );


    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);
    Info << "correcting Nut" << endl;
    this->correctNut(S2, F2);
    Info << "correcting nonlinear stress" << endl;

    correctNonlinearStress(gradU);

    tgradU.clear();
    Info << this->fK_ << endl;
    Info << this->fOmega_ << endl;
    Info << "End" << endl;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
