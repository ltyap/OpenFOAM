/*---------------------------------------------------------------------------*\
Date: June 24, 2005
Author: Eugene de Villiers
Source: http://www.cfd-online.com/Forums/openfoam-solving/58043-les-2.html#post187619
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    perturbU

Description
    initialise channel velocity with superimposed streamwise streaks.
    To be used to force channelOodles to transition and reach a fully
    developed flow sooner.

    Reads in perturbUDict.

    EdV from paper:
        Schoppa, W. and Hussain, F.
        "Coherent structure dynamics in near wall turbulence",
        Fluid Dynamics Research, Vol 26, pp119-139, 2000.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "wallDist.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // These values should be edited according to the case
    // Works for annular pipe flows with origin at the center
    // And flow is in the x-direction
    // h is the half-height of the annular gap
    // Ubar should be set in transportProperties

    IOdictionary perturbDict
    (
        IOobject
        (
            "perturbUAnnularDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const scalar RetauInner(readScalar(perturbDict.lookup("RetauInner")));
    const scalar RetauOuter(readScalar(perturbDict.lookup("RetauOuter")));
    const scalar eta(readScalar(perturbDict.lookup("eta")));
    const scalar Ro(readScalar(perturbDict.lookup("Ro")));
    const bool setBulk(readBool(perturbDict.lookup("setBulk")));
    const bool perturb(readBool(perturbDict.lookup("perturb")));

    // Perturbation properties
    const scalar duPlusCoeff{readScalar(perturbDict.lookup("duPlusCoeff"))};
    scalar betaPlus{readScalar(perturbDict.lookup("betaPlus"))};
    scalar alphaPlus{readScalar(perturbDict.lookup("alphaPlus"))};
    const scalar sigma{readScalar(perturbDict.lookup("sigma"))};
    const scalar epsilonCoeff{readScalar(perturbDict.lookup("epsilonCoeff"))};

    const scalar Ri = eta*Ro;

    Info<< "Radius ratio       = " << eta << nl
        << "Outer radius       = " << Ro << nl
        << "Inner radius       = " << Ri << nl
        << "Inner Re(tau) = " << RetauInner << nl
        << "Outer Re(tau) = " << RetauOuter << nl
        << "Set bulk flow             = " << Switch(setBulk) << nl
        << "Perturb flow              = " << Switch(perturb) << nl
        << endl;


    if (!setBulk && !perturb)
    {
        FatalErrorIn(args.executable())
            << "At least one of setBulk or perturb needs to be set"
            << " to do anything to the velocity"
            << exit(FatalError);
    }

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );
    Info<< "Reading U" << endl;

    volVectorField U(Uheader, mesh);

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
    	"nu",
    	dimViscosity,
        transportProperties.lookup("nu")
    );
    dimensionedVector Ubar
    (
    	"Ubar",
    	dimVelocity,
        transportProperties.lookup("Ubar")
    );


    Info<< "nu      = " << nu << endl;
    Info<< "Ubar    = " << Ubar << endl;
    const scalar delta = 0.5*(Ro-Ri);
    const scalar utauInner = RetauInner*nu.value()/delta;
    const scalar utauOuter = RetauOuter*nu.value()/delta;
    Info<< "Annular half width gap (delta) = " << delta << endl;
    Info<< "Inner u(tau)  = " << utauInner << endl;
    Info<< "Outer u(tau)  = " << utauOuter << endl;

    //spanwise wavenumber: spacing z+ = 200
    betaPlus = 2.0*constant::mathematical::pi*(1.0/betaPlus);
    // const scalar sigma = 0.00055;
    //streamwise wave number: spacing x+ = 500
    alphaPlus = 2.0*constant::mathematical::pi*(1.0/alphaPlus);
    const scalar epsilon = epsilonCoeff*Ubar.value()[0];

    // Random number generator
    Random perturbation(1234567);

    const vectorField& centres = mesh.C();

    const vector zDir = Ubar.value()/Foam::mag(Ubar.value()); // streamwise vector
    Info<< "zDir = " << zDir << endl;
    wallDist yWD(mesh);
    const volVectorField& rDir = yWD.n(); //wall normal vectors

    const scalar RobyRi = Ro/Ri;
    const scalar diffSqrRoRi = Foam::sqr(Ro) - Foam::sqr(Ri);
    const scalar coeff1 = 2 * (diffSqrRoRi)
        / (Foam::pow4(Ro) - Foam::pow4(Ri) - (Foam::sqr(diffSqrRoRi)/(Foam::log(RobyRi))));


    forAll(centres, celli)
    {
        // add a small (+/-20%) random component to enhance symetry breaking
        scalar deviation=1.0 + 0.2*perturbation.scalarNormal();

        const vector& cCentre = centres[celli];


        const scalar r = Foam::mag(cCentre & rDir[celli]); // distance from origin
        // const scalar r = ::sqrt(::sqr(cCentre.y()) + ::sqr(cCentre.z()));

        const scalar theta = Foam::atan(cCentre.y()/cCentre.z());

        const scalar yWall = min(r-Ri, Ro-r); 

        const scalar Retau = r < (delta + Ri) ? RetauInner : RetauOuter;
        const scalar utau = r < (delta + Ri) ? utauInner : utauOuter;

        vector thetaDir = r < (delta + Ri) ? rDir[celli]^zDir : zDir^rDir[celli];
        thetaDir /= mag(thetaDir);
/*
        vector tangential
        (
            0, cCentre.y(), cCentre.z()
        );
        tangential = tangential ^ vector(1,0,0);
        tangential = tangential/mag(tangential);

        Info<< "thetaDir" << thetaDir <<endl;
        Info<< "thetaDir alternate" << tangential << endl;
*/

        //wall normal circulation
        const scalar duplus = Ubar.value()[0]*duPlusCoeff/utau;

        const scalar yplus = yWall*Retau/delta; 

        const scalar zplus = (cCentre & zDir)*Retau/delta;
        const scalar rplus = r*Retau/delta; 
        // Info << "rplus = " << rplus << endl;

        if (setBulk)
        {
            // laminar parabolic profile
            U[celli] = vector::zero;

            const scalar rByRi = r/Ri;
            U[celli] =
                coeff1*Ubar.value()
                * (Foam::sqr(Ri) - Foam::sqr(r) + diffSqrRoRi * Foam::log(rByRi)/Foam::log(RobyRi));
        }

        if (perturb)
        {
            // streak streamwise velocity
            U[celli] +=
                zDir*(utau * duplus/2.0) * (yplus/30.0)
                * Foam::exp(-sigma * Foam::sqr(yplus) + 0.5)
                * Foam::cos(betaPlus*rplus*theta)*deviation;

            // streak spanwise perturbation
            U[celli] +=
               thetaDir * epsilon
              * Foam::sin(alphaPlus*zplus)
              * yplus
              * Foam::exp(-sigma*Foam::sqr(yplus))
              * deviation;
         }
        
    }

    Info<< "Writing modified U field to " << runTime.timeName() << endl;
    U.write();

    Info<< endl;
    

    return(0);
}


// ************************************************************************* //
