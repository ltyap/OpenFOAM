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
    Foam::RASModels::QCR24kOSST_PANS

Description
    Specialisation for RAS of the generic kOmegaSSTBase base class.
    Implements the QCR model of Tamaki, Y. and Kawai, S (QCR2024)
    For more information, see Description of kOmegaSSTBase.H

See also
    Foam::kOmegaSST

SourceFiles
    QCR24kOSST_PANS.C

\*---------------------------------------------------------------------------*/

#ifndef QCR24kOSST_H
#define QCR24kOSST_H

#include "kOmegaSSTPANSBase.H"
#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSSTQCR Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class QCR24kOSST_PANS
:
    public Foam::kOmegaSSTPANS
    <
        nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
{

protected:
    virtual void correctNonlinearStress(const volTensorField& gradU);

    dimensionedScalar Cr1_;
    dimensionedScalar Cr2_;




public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("QCR24kOSST_PANS");


    // Constructors

        //- Construct from components
        QCR24kOSST_PANS
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

    virtual bool read();

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();

    //- Destructor
    virtual ~QCR24kOSST_PANS()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "QCR24kOSST_PANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
