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

Class
    Foam::RASModels::kOmegaSST_HTPANS

Description
    Partially-Average Navier-Stokes model based on the k-omega-SST RAS model.
    Details of the coefficients can be found in the kOmegaSSTBase class.

    The default model coefficients are
    \verbatim
        kOmegaSSTCoeffs
        {
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;
        }
    \endverbatim
    
    15/08/2025:
    Note that the Prandtl number alphaOmega2 in the blending function F1
    is not recalibrated by fOmega/fK as described in
    \verbatim
        F.S. Pereira, G. Vaz, L. Eça, S.S. Girimaji,
        Simulation of the flow around a circular cylinder at Re=3900 with Partially-Averaged Navier-Stokes equations,
        International Journal of Heat and Fluid Flow,
        Volume 69,
        2018,
        Pages 234-246
    \endverbatim

See also
    Foam::kOmegaSST

SourceFiles
    kOmegaSST_HTPANS.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTHTPANS_H
#define kOmegaSSTHTPANS_H

#include "kOmegaSSTPANS.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSSTPANS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaSST_HTPANS
:
    public kOmegaSSTPANS<BasicMomentumTransportModel>
{
protected:
    // Protected data
            IOdictionary transportProperties_;

        // Model coefficients for HT model

            dimensionedScalar Clambda_;
            dimensionedScalar Cp1_;
            dimensionedScalar Cp2_;
            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;
            dimensionedScalar Cgamma_;
            dimensionedScalar Ceps2_;

            dimensionedScalar sigmakTheta_;
            dimensionedScalar sigmaEpsTheta_;
  
            dimensionedScalar Pr_;

            //- Lower limit of kTheta
            dimensionedScalar kThetaMin_;

            //- Lower limit of epsilon
            dimensionedScalar epsilonThetaMin_;


            volScalarField kTheta_;
            volScalarField epsilonTheta_;

            // Thermal diffusivity
            volScalarField Dt_;
            volScalarField D_;

            volScalarField fKTheta_;
            volScalarField fEpsTheta_;

            // Functions for HT model
            tmp<volScalarField> yStar() const;
            tmp<volScalarField> Rt() const;
            tmp<volScalarField> fEps2() const;
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
    TypeName("kOmegaSST_HTPANS");


    // Constructors

        //- Construct from components
        kOmegaSST_HTPANS
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        kOmegaSST_HTPANS(const kOmegaSST_HTPANS&) = delete;


    //- Destructor
    virtual ~kOmegaSST_HTPANS()
    {}

    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

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

        virtual tmp<volScalarField> fKTheta() const
        {
            return fKTheta_;
        }

        virtual tmp<volScalarField> fEpsTheta() const
        {
            return fEpsTheta_;
        }


        //- Solve the thermal equations and correct the turbulence thermal diffusivity
        virtual void correctThermal();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

    // Member Operators

        //- No copy assignment
        void operator=(const kOmegaSST_HTPANS&) = delete;




};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSST_HTPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
