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
    Foam::RASModels::AKNPANS

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
            Cmu         0.09;
            Ceps1       1.5;
            Ceps2       1.9;
            sigmaEps    1.4;
            sigmaK      1.4;
        }
    \endverbatim

SourceFiles
    AKNPANS.C

\*---------------------------------------------------------------------------*/

#ifndef AKNPANS_H
#define AKNPANS_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class AKNPANS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class AKNPANS
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Ceps3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            volScalarField fK_;
            volScalarField fEps_;


    // Protected Member Functions
        tmp<volScalarField> yStar() const;
        tmp<volScalarField> Rt() const;
        tmp<volScalarField> fMu() const;
        tmp<volScalarField> fEps2() const;
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("AKNPANS");


    // Constructors

        //- Construct from components
        AKNPANS
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
        AKNPANS(const AKNPANS&) = delete;

    //- Destructor
    virtual ~AKNPANS()
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

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

    // Member Operators

        //- No copy assignment
        void operator=(const AKNPANS&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AKNPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
