/*---------------------------------------------------------------------------*\
Application
    postLES

Description
	Calculates second moment statistics based on UMean, together with production and dissipation
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"

	#include "readTransportProperties.H"
	
    #include "readFields.H"   
	
	volScalarField Tp2Mean
	(
		IOobject
        (   
            "Tp2Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		TTMean - (TMean * TMean)
	);
    volVectorField UPrimeTMean
	(
		IOobject
        (   
            "UPrimeTMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		uiTMean - (UMean * TMean)
	);
	// volScalarField production
	// (
	// 	IOobject
    //     (   
    //         "production",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
	// 	-(Up2Mean && djuiMean)
	// );
	
	volScalarField scalarDiss
	(
		IOobject
        (   
            "scalarDiss",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		-alpha*(djTdjTMean - (djTMean & djTMean))
	);

	
    Info<< "Writing Tp2Mean to " << runTime.timeName() << endl;
 	Tp2Mean.write();

    Info<< "Writing UPrimeTMean to " << runTime.timeName() << endl;
 	UPrimeTMean.write();

    Info<< "Writing scalarDiss to " << runTime.timeName() << endl;
 	scalarDiss.write();



    Info<< endl;
    

    return(0);
}


// ************************************************************************* //
