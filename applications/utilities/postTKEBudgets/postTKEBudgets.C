/*---------------------------------------------------------------------------*\
Application
    postTKEBudgets

Description
	Calculates second moment statistics based on UMean, together with production and dissipation
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculate full set of TKE budget terms"
    );

    argList::noParallel();
    timeSelector::addOptions();

	#include "setRootCase.H"
	#include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

	#include "createMesh.H"
	#include "readTransportProperties.H"
	
    forAll(timeDirs, timeI)
    {
        #include "readFields.H"   
            
        volSymmTensorField Up2Mean
        (
            IOobject
            (   
                "Up2Mean",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            uiujMean - symm(UMean * UMean)
        );
        
        volScalarField production
        (
            IOobject
            (   
                "production",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -(Up2Mean && djuiMean)
        );
        
        volScalarField dissipation
        (
            IOobject
            (   
                "dissipation",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -nu*(djuidjuiMean - (djuiMean && djuiMean))
        );

        volScalarField diff
        (
            IOobject
            (   
                "diff",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            nu*(uiLuiMean - (UMean & LuiMean))
        );

        volScalarField pDiff
        (
            IOobject
            (   
                "pDiff",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -(uiGcpMean - (UMean & GcpMean))
        );

        volScalarField conv
        (
            IOobject
            (   
                "conv",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -(uiCuiMean - (UMean & CuiMean))
        );

        Info<< "Writing Up2Mean to " << runTime.timeName() << endl;
        Up2Mean.write();
        Info<< "Writing production to " << runTime.timeName() << endl;
        production.write();
        Info<< "Writing dissipation to " << runTime.timeName() << endl;
        dissipation.write();

        Info<< "Writing diff to " << runTime.timeName() << endl;
        diff.write();
        Info<< "Writing pDiff to " << runTime.timeName() << endl;
        pDiff.write();
        Info<< "Writing conv to " << runTime.timeName() << endl;
        conv.write();


        Info<< endl;
    }
 
    

    return(0);
}


// ************************************************************************* //
