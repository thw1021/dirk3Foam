    Info<< "Reading transportProperties\n" << endl;

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
        transportProperties
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;


    

/*---------------------------------------------------------------*/
Info<< "Reading/calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);
/*---------------------------------------------------------------*/

// R1 residual vector field of the 1st stage
volVectorField R1
(
    IOobject( "R1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mesh,
    dimensionedVector("R1", dimAcceleration, vector::zero)
);

// R2 residual vector field of the 2nd stage
volVectorField R2
(
    IOobject( "R2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mesh,
    dimensionedVector("R2", dimAcceleration, vector::zero)
);

// R1 residual vector field of the 1st stage
surfaceVectorField Rf1
(
    IOobject( "Rf1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mesh,
    dimensionedVector("Rf1", dimAcceleration, vector::zero)
);

// R2 residual vector field of the 2nd stage
surfaceVectorField Rf2
(
    IOobject( "Rf2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ),
    mesh,
    dimensionedVector("Rf2", dimAcceleration, vector::zero)
);

volVectorField Uold = U;

