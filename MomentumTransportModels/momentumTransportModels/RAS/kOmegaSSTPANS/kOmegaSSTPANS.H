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
    Foam::RASModels::kOmegaSSTPANS

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
    kOmegaSSTPANS.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTPANS_H
#define kOmegaSSTPANS_H

#include "kOmegaSSTPANSBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSSTPANS Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaSSTPANS
:
    public Foam::kOmegaSSTPANS
    <
        eddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
{

  
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTPANS");


    // Constructors

        //- Construct from components
        kOmegaSSTPANS
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
        kOmegaSSTPANS(const kOmegaSSTPANS&) = delete;


    //- Destructor
    virtual ~kOmegaSSTPANS()
    {}

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
