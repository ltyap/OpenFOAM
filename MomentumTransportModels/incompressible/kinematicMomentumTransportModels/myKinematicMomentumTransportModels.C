/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "kinematicMomentumTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "AbeKondohNaganoKE.H"
makeRASModel(AbeKondohNaganoKE);

#include "AKNHT.H"
makeRASModel(AKNHT);

#include "AKN_DWX.H"
makeRASModel(AKN_DWX);
/*
#include "AbeKondohNaganoKE_HT.H"
makeRASModel(AbeKondohNaganoKE_HT);

#include "AbeKondohNaganoKE_HTPANS.H"
makeRASModel(AbeKondohNaganoKE_HTPANS);

#include "AKNPANS.H"
makeRASModel(AKNPANS);

#include "AKN_MMPANS.H"
makeRASModel(AKN_MMPANS);

#include "kOmegaSST_HTPANS.H"
makeRASModel(kOmegaSST_HTPANS);

*/
#include "kOmegaSSTPANS.H"
makeRASModel(kOmegaSSTPANS);

#include "kOmegaPANS.H"
makeRASModel(kOmegaPANS);

/*
#include "QCRkOSST.H"
makeRASModel(QCRkOSST);

#include "QCR24kOSST_PANS.H"
makeRASModel(QCR24kOSST_PANS);

#include "QCRkOSSTPANS.H"
makeRASModel(QCRkOSSTPANS);
*/

// -------------------------------------------------------------------------- //
// PANS models
// -------------------------------------------------------------------------- //

// ************************************************************************* //
