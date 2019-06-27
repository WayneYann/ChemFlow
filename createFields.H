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
 volScalarField Ydefault
 (
     IOobject
     (
         "Ydefault",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE,
         false
     ),
     mesh,
     dimensionedScalar("Ydefault", dimless, 1)
 );

 Ydefault.write();

 volScalarField p
 (
     IOobject
     (
         "p",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE,
         false
     ),
     mesh,
     dimensionedScalar("p", dimPressure, p0)
 );

 p.write();

 volScalarField T
 (
     IOobject
     (
         "T",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE,
         false
     ),
     mesh,
     dimensionedScalar("T", dimTemperature, 300.0)
 );

 T.write();
}