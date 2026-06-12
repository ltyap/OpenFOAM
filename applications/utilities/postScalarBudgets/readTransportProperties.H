Info << "Reading transportProperties" << endl;

IOdictionary transportProperties
( 
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar Pr
(
    "Pr",
    dimless,
    transportProperties.lookup("Pr")
);

dimensionedScalar alpha(nu/Pr);