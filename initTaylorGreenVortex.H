
// This file is only used to test dirk3Foam solver

// *** This initialization is used to test Taylor-Green vortex case *********************************
    Info<< "Calculating initial conditions" << endl;
    // Initial conditions based on analytical solution at t=0
    forAll(U,cellI)
    {
        U[cellI].x()=-Foam::cos(U.mesh().C()[cellI].x()) * Foam::sin(U.mesh().C()[cellI].y());
        U[cellI].y()= Foam::sin(U.mesh().C()[cellI].x()) * Foam::cos(U.mesh().C()[cellI].y());
        p[cellI]=-0.25*(Foam::cos(2*U.mesh().C()[cellI].x()) + Foam::cos(2*U.mesh().C()[cellI].y()));
    }
// ***************************************************************************************************