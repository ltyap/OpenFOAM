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
    for incompressible flows, with heat transfer (scalar transport)
	
    Wall boundary conditions for the momentum field are:

        k       = 1e-10
        epsilon = epsilonWallFunction
    
    Wall boundary conditions for the thermal field are:
        
        kTheta          = 1e-10
        epsilonTheta    = epsilonTheta2WallFunction

	
	References:

    Abe K., Kondoh T., Nagano Y., 1995. A new turbulence model for predicting
	fluid flow and heat transfer in separating and reattaching flows-II. 
    Thermal field calculations. Int. J. Heat and Mass Transfer 38(8), pp 1467-1481



    The default model coefficients for the thermal model are given as:
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

See also
    Foam::RAS::AbeKondohNaganoKE

SourceFiles
    AKNHT.C

\*---------------------------------------------------------------------------*/

#ifndef AKNHT_H
#define AKNHT_H

#include "AbeKondohNaganoKE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class AKN Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class AKNHT
:
    public AbeKondohNaganoKE<BasicMomentumTransportModel>
{

protected:

    // Protected data

        // Model coefficients
            IOdictionary transportProperties_;

            dimensionedScalar Clambda_;
            dimensionedScalar Cp1_;
            dimensionedScalar Cp2_;
            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;

            dimensionedScalar Cgamma_;

            dimensionedScalar Pr_;

            dimensionedScalar sigmakTheta_;
            dimensionedScalar sigmaEpsTheta_;

            //- Lower limit of kTheta
            dimensionedScalar kThetaMin_;

            //- Lower limit of epsilon
            dimensionedScalar epsilonThetaMin_;


        // Fields

            volScalarField kTheta_;
            volScalarField epsilonTheta_;

            volScalarField Dt_;
            volScalarField D_;

        // Protected Member Functions
            tmp<volScalarField> fp1() const;
            tmp<volScalarField> fd1() const;

            tmp<volScalarField> fd2() const;
            tmp<volScalarField> R() const;
            tmp<volScalarField> fAlpha() const;
            
            virtual void correctAlphat();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("AKNHT");


    // Constructors

        //- Construct from components
        AKNHT
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
        AKNHT(const AKNHT&) = delete;

    //- Destructor
        virtual ~AKNHT()
        {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        tmp<volScalarField> D()
        {
            return D_;
        }

        tmp<volScalarField> Dt()
        {
            return Dt_;
        }

        //- Return the effective diffusivity for kTheta
        tmp<volScalarField> DkThetaEff() const
        {
            return volScalarField::New
            (
                "DkThetaEff",
                (Dt_/sigmakTheta_) + D_
            );
        }
        
        //- Return the effective diffusivity for epsilonTheta
        tmp<volScalarField> DepsilonThetaEff() const
        {
            return volScalarField::New
            (
                "DepsilonThetaEff",
                (Dt_/sigmaEpsTheta_) + D_
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

        //- Solve the thermal equations and correct the turbulence thermal diffusivity
        virtual void correctThermal();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

    // Member Operators

        //- No copy assignment
        void operator=(const AKNHT&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AKNHT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
