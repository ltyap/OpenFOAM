IOobject UMeanheader
(
    "UMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading UMean" << endl;

volVectorField UMean(UMeanheader, mesh);

IOobject uiujMeanheader
(
    "uiujMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading uiujMean" << endl;

volSymmTensorField uiujMean(uiujMeanheader, mesh);

IOobject djuidjuiMeanheader
(
    "djuidjuiMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading djuidjuiMean" << endl;

volScalarField djuidjuiMean(djuidjuiMeanheader, mesh);

IOobject djuiMeanheader
(
    "djuiMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading djuiMean" << endl;

volTensorField djuiMean(djuiMeanheader, mesh);

volVectorField CuiMean
(
    IOobject
    (
        "CuiMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);

volScalarField uiCuiMean
(
    IOobject
    (
        "uiCuiMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);

volVectorField LuiMean
(
    IOobject
    (
        "LuiMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);

volScalarField uiLuiMean
(
    IOobject
    (
        "uiLuiMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);

volVectorField GcpMean
(
    IOobject
    (
        "GcpMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);

volScalarField uiGcpMean
(
    IOobject
    (
        "uiGcpMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);















