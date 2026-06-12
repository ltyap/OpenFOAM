IOobject UMeanheader
(
    "UMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading UMean" << endl;

volVectorField UMean(UMeanheader, mesh);

IOobject TMeanheader
(
    "TMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading TMean" << endl;

volScalarField TMean(TMeanheader, mesh);

IOobject TTMeanheader
(
    "TTMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading TTMean" << endl;

volScalarField TTMean(TTMeanheader, mesh);

IOobject djTdjTMeanheader
(
    "djTdjTMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading djTdjTMean" << endl;

volScalarField djTdjTMean(djTdjTMeanheader, mesh);

IOobject djTMeanheader
(
    "djTMean",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ
);
Info<< "Reading djTMean" << endl;

volVectorField djTMean(djTMeanheader, mesh);

volVectorField uiTMean
(
    IOobject
    (
        "uiTMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    ),
    mesh
);
