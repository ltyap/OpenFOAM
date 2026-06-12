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
    Foam::RASModels::QCRkOSSTPANS

Description
    Specialisation for PANs of the generic kOmegaSSTPANSBase base class.
    Implements the QCR model of Spalart (QCR2000)
    For more information, see Description of kOmegaSSTPANSBase.H

See also
    Foam::kOmegaSSTPANS

SourceFiles
    QCRkOSSTPANS.C

\*---------------------------------------------------------------------------*/

#ifndef QCRkOSSTPANS_H
#define QCRkOSSTPANS_H

#include "kOmegaSSTPANSBase.H"
#include "RASModel.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSSTPANSQCR Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class QCRkOSSTPANS
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

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("QCRkOSSTPANS");


    // Constructors

        //- Construct from components
        QCRkOSSTPANS
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );
    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

    //- Destructor
    virtual ~QCRkOSSTPANS()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "QCRkOSSTPANS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
