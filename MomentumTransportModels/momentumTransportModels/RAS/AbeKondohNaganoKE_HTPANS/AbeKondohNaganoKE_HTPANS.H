/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::RASModels::AKNHT

Group
    grpRASTurbulence

Description
    Abe, Kondoh and Nagano's low-Reynolds number k-epsilon turbulence model
    for incompressible flows
	
    Wall boundary conditions are:

        k       = 1e-10
        epsilon = epsilonWallFunction
	
	
	References:

    Abe K., Kondoh T., Nagano Y., 1994. A new turbulence model for predicting
	fluid flow and heat transfer in separating and reattaching flows-I. Flow
	field calculations. Int. J. Heat and Mass Transfer 37(1), pp 138-151

    The default model coefficients are given as:
    \verbatim
        AbeKondohNaganoKECoeffs
        {
            Clambda              0.10;
            Cp1                  1.9;
            Cp2                  0.6;
            Cd1                  2.0;
            Cd2                  0.9;
            sigmaEpsTheta        1.6;
            sigmaKTheta          1.6;
        }
    \endverbatim

SourceFiles
    AKNHT.C

\*---------------------------------------------------------------------------*/

#ifndef AbeKondohNaganoKE_HTPANS_H
#define AbeKondohNaganoKE_HTPANS_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class AbeKondohNaganoKE_HTPANS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class AbeKondohNaganoKE_HTPANS
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{

protected:

    // Protected data
            IOdictionary transportProperties_;

        // Model coefficients for turbulence model

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Ceps3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Model coefficients for HT model

            dimensionedScalar Clambda_;
            dimensionedScalar Cp1_;
            dimensionedScalar Cp2_;
            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;
            dimensionedScalar Cgamma_;
            dimensionedScalar sigmakTheta_;
            dimensionedScalar sigmaEpsTheta_;

            dimensionedScalar Pr_;

            //- Lower limit of kTheta
            dimensionedScalar kThetaMin_;

            //- Lower limit of epsilon
            dimensionedScalar epsilonThetaMin_;


        // Fields
            volScalarField k_;
            volScalarField epsilon_;
            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            volScalarField kTheta_;
            volScalarField epsilonTheta_;

            // Thermal diffusivity
            volScalarField Dt_;
            volScalarField D_;

            volScalarField fK_;
            volScalarField fEps_;

            volScalarField fKTheta_;
            volScalarField fEpsTheta_;

        // Protected Member Functions

            // Functions for turbulence model
            tmp<volScalarField> yStar() const;
            tmp<volScalarField> Rt() const;
            tmp<volScalarField> fMu() const;
            tmp<volScalarField> fEps2() const;
            virtual void correctNut();

            // Functions for HT model
            tmp<volScalarField> fp1() const;
            tmp<volScalarField> fd1() const;

            tmp<volScalarField> fd2() const;
            tmp<volScalarField> R() const;
            tmp<volScalarField> fAlpha() const;
            
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("AbeKondohNaganoKE_HTPANS");


    // Constructors

        //- Construct from components
        AbeKondohNaganoKE_HTPANS
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- No copy construct
        AbeKondohNaganoKE_HTPANS(const AbeKondohNaganoKE_HTPANS&) = delete;

    //- Destructor
        virtual ~AbeKondohNaganoKE_HTPANS()
        {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (
                    (this->nut_/(sqr(fK_)*sigmak_/fEps_)) 
                    + this->nu()
                )
            );
        }
        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (
                    (this->nut_/(sqr(fK_)*sigmaEps_/fEps_))
                    + this->nu()
                )
            );
        }
       
        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        virtual tmp<volScalarField> D() const
        {
            return D_;
        }

        virtual tmp<volScalarField> Dt() const
        {
            return Dt_;
        }

        //- Return the effective diffusivity for kTheta
        tmp<volScalarField> DkThetaEff() const
        {
            return volScalarField::New
            (
                "DkThetaEff",
                (Dt_/(sigmakTheta_*sqr(fKTheta_)/fEpsTheta_)) + D_
            );
        }
        
        //- Return the effective diffusivity for epsilonTheta
        tmp<volScalarField> DepsilonThetaEff() const
        {
            return volScalarField::New
            (
                "DepsilonThetaEff",
                (Dt_/(sigmaEpsTheta_*sqr(fKTheta_)/fEpsTheta_)) + D_
            );
        }


        //- Return the turbulence thermal variance
        virtual tmp<volScalarField> kTheta() const
        {
            return kTheta_;
        }

        //- Return the turbulence thermal variance dissipation rate
        virtual tmp<volScalarField> epsilonTheta() const
        {
            return epsilonTheta_;
        }

        //- Solve the four turbulence equations and correct  both the turbulence viscosity and turbulent diffusivity
        virtual void correct();

    // Member Operators

        //- No copy assignment
        void operator=(const AbeKondohNaganoKE_HTPANS&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AbeKondohNaganoKE_HTPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
