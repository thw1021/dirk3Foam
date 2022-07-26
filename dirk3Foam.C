/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    dirk3Foam


Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    DIRKRK3-3: Diagonally Implict Runge-Kutta (3rd order, 3 stages) Time Scheme

    Ref 1: Valerio D’Aless, Lorenzo Binci, Sergio Montelpare, Renato Ricci.
            "On the development of OpenFOAM solvers based on explicit and implicit 
             high-order Runge–Kutta schemes for incompressible flows with heat transfer",
            Computer Physics Communications, 2018(222):14-30.

    Ref 2: E.M.J. Komen, E.M.A. Frederix, T.H.J. Coppen, V. D’Alessandro, J.G.M. Kuerten.
            "Analysis of the numerical dissipation rate of different Runge–Kutta and
            velocity interpolation methods in an unstructured collocated finite volume method in OpenFOAM.""
            Computer Physics Communications 253 (2020) 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createControl.H"      
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        Info<< "SDIRK 3--3 "  << nl << endl;

    //#include "readCoeff.H"
    #include "sandersCoeff.H"

    Info<< "Sanderse RK -- Scheme 1 -- Coeff. a11 "<< a11  << nl << endl;
    Info<< "Sanderse RK -- Scheme 1 -- Coeff. a21 "<< a21  << nl << endl;
    Info<< "Sanderse RK -- Scheme 1 -- Coeff. a22 "<< a22  << nl << endl;
    Info<< "Sanderse RK -- Scheme 1 -- Coeff. a31 "<< a31  << nl << endl;
    Info<< "Sanderse RK -- Scheme 1 -- Coeff. a32 "<< a32  << nl << endl;
    Info<< "Sanderse RK -- Scheme 1 -- Coeff. a33 "<< a33  << nl << endl;


    dimensionedVector zeroU ( "zeroU", dimensionSet(0, 1, -2, 0, 0, 0, 0),vector(0, 0, 0) );

    Info<< "\nStarting time loop\n" << endl;


    // Note: only include the file for testing Taylor-Green vortex case
    #include "initTaylorGreenVortex.H"
    //

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
 

        #include "CourantNo.H"
        
//----------------------------------------------------------------------------
       dimensionedScalar dt = runTime.deltaT();
       Uold = U;    

    // ***********  stage 1 **************
        Info << nl << "Stage 1 "  << nl << endl;

        // * 1) solve disretized momentum equation
        fvVectorMatrix UEqn(fvm::div(phi, U) - fvm::laplacian(nu, U) + fvm::Sp(1.0 / (a11*dt), U)); // Eqn. (25)
        solve(UEqn == Uold / (a11 * dt) - fvc::grad(p));                                              // Eqn. (19)

        while (piso.correct())
        {
        // * 2) calculate mass fluxes
            volScalarField rAU("rAU", 1.0 / UEqn.A());  // Please note, here 1/ap has been updated to be (1 / ap_hat)
            volVectorField HbyA(constrainHbyA(rAU * ( Uold / (dt * a11) + UEqn.H() ), U, p));          // Eqn. 21 (modified from pisoFoam)
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

            adjustPhi(phiHbyA, U, p);
            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // * 3) solve pressure equation
            fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));
            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

            if (piso.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
                // phi -= pEqn.flux(); // Seems the same as the above                
            }

            // * 4) perform velocity correction
            //U = rAU * (Uold / (dt*a11)  + UEqn.H() - fvc::grad(p));  // Eqn. 22  The same as the line below
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
     
        }

        // The tricky part: retrieve ap form ap_hat            
        R1 = -(UEqn.A()-(1.0/(a11*dt)))*U + UEqn.H() - fvc::grad(p);            // Eqn. 17
        // face residual vector field  -- Eqn. 23
        Rf1 = ( (fvc::interpolate(U) - fvc::interpolate(Uold)) / (a11 * dt) );    // Eqn. (23)

        // ***********  stage 2 **************
        Info<< nl << "Stage 2 "  << nl << endl;

        // * 1) solve disretized momentum equation
        //fvVectorMatrix UEqn( fvm::div(phi, U) - fvm::laplacian(nu, U) + fvm::Sp(1.0 / (a22*dt), U) );// Eqn. (25) + (19)
        solve(UEqn == Uold / (a22 * dt) - fvc::grad(p) + (a21 * R1) / a22 );


        while (piso.correct())
        {
        // * 2) calculate mass fluxes
            volScalarField rAU("rAU", 1.0/UEqn.A());          
            volVectorField HbyA(constrainHbyA(rAU * (Uold / (dt * a22) + UEqn.H() ), U, p));         // Eqn. (21)
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            phiHbyA = phiHbyA + ((fvc::interpolate(rAU) * Rf1 * a21 / a22) & mesh.Sf());            // Eqn. (21)
            HbyA += rAU * (a21 * R1 / a22);

            adjustPhi(phiHbyA, U, p);
            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

        // * 3) solve pressure equation
            fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));                       // En. (20)
            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

            if (piso.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
                //phi -= pEqn.flux(); The same as the above
            }          

        // * 4) perform velocity correction
            //U = rAU * (Uold / (a22 * dt) + UEqn.H() - fvc::grad(p) + a21 / a22 * R1);             // Eqn. (22)
            U = HbyA - rAU * fvc::grad(p);
            U.correctBoundaryConditions();
  
        }

        // * 5) calculate the residual for next stage
        // The tricky part: retrieve ap form ap_hat            
        R2  = -(UEqn.A() - 1.0/(a22*dt)) * U  + UEqn.H() - fvc::grad(p);                            // Eqn (17)
        // face residual vector field  --- Eqn (23)
        Rf2 = (fvc::interpolate(U) - fvc::interpolate(Uold)) / (a22 * dt) - (a21 * Rf1) / a22 ;    // Eq. (23)

        // ***********  stage 3 **************
        Info<< nl << "Stage 3 "  << nl << endl;

        // * 1) solve disretized momentum equation
        solve(UEqn == Uold / (a33 * dt) - fvc::grad(p) + (a31 * R1 + a32 * R2) / a33); // Eqn. (25) + (19)
        //U.correctBoundaryConditions();
        while (piso.correct())
        {
        // * 2) calculate mass fluxes
            volScalarField rAU("rAU", 1.0/UEqn.A());          
            volVectorField HbyA(constrainHbyA(rAU * (Uold / (dt * a33) + UEqn.H()), U, p));
            surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
            phiHbyA = phiHbyA + ((fvc::interpolate(rAU) * (Rf1 * a31 + Rf2 * a32) / a33) & mesh.Sf()); // Eqn. 21
            HbyA += rAU * (a31 * R1 + a32 * R2) / a33;

            adjustPhi(phiHbyA, U, p);
            //Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

        // * 3) solve pressure equation

            fvScalarMatrix pEqn(fvm::laplacian(rAU, p) == fvc::div(phiHbyA));
            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
            if (piso.finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
                //phi -= pEqn.flux(); // the same as the above
            }

        // * 4) perform velocity correction
            #include "continuityErrs.H"
            //U = rAU * (Uold / (a33 * dt) + UEqn.H() - fvc::grad(p) + a31 / a33 * R1 + a32 / a33 * R2);
            // This line is the same as the above line but saves some flops.
            U = HbyA - rAU * fvc::grad(p);
            U.correctBoundaryConditions();

        }

        phi = ( fvc::interpolate(U)  & mesh.Sf());

        //----------------------------------------------------------------------------

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
