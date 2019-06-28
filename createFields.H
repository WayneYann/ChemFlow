// Create a single cell mesh to hold OpenFOAM thermo
// Note that the actual geometry does not affect anything
Foam::argList args(argc, argv);
Foam::Time runTime(Foam::Time::controlDictName, args);
Foam::labelList owner(6, Foam::label(0));
Foam::labelList neighbour(0);

Foam::pointField points(8);
points[0] = Foam::vector(0, 0, 0);
points[1] = Foam::vector(1, 0, 0);
points[2] = Foam::vector(1, 1, 0);
points[3] = Foam::vector(0, 1, 0);
points[4] = Foam::vector(0, 0, 1);
points[5] = Foam::vector(1, 0, 1);
points[6] = Foam::vector(1, 1, 1);
points[7] = Foam::vector(0, 1, 1);

const Foam::cellModel& hexa = *(Foam::cellModeller::lookup("hex"));
Foam::faceList faces = hexa.modelFaces();

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::READ_IF_PRESENT
    ),
    Foam::xferMove<Foam::Field<Foam::vector>>(points),
    faces.xfer(),
    owner.xfer(),
    neighbour.xfer()
);

Foam::List<Foam::polyPatch*> patches(1);

patches[0] = new Foam::emptyPolyPatch
(
    "boundary",
    6,
    0,
    0,
    mesh.boundaryMesh(),
    Foam::emptyPolyPatch::typeName
);

mesh.addFvPatches(patches);

{
    Foam::volScalarField Ydefault
    (
        Foam::IOobject
        (
            "Ydefault",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh,
        Foam::dimensionedScalar("Ydefault", dimless, 1)
    );

    Ydefault.write();

    Foam::volScalarField p
    (
        Foam::IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh,
        Foam::dimensionedScalar("p", dimPressure, 101325)
    );

    p.write();

    Foam::volScalarField T
    (
        Foam::IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE,
            false
        ),
        mesh,
        Foam::dimensionedScalar("T", dimTemperature, 300)
    );

    T.write();
}