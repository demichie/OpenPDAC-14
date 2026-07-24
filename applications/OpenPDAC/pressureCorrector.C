/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenPDAC.
    This file was derived from the multiphaseEuler solver in OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "OpenPDAC.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::pressureCorrector()
{

    const Switch logPressureExtrema(
        pimple.dict().lookupOrDefault<Switch>("logPressureExtrema", false));

    if (forceFinalPimpleIter_ && !pimple.finalIter())
    {
        Info << "Bypass p_rgh residual history size before append = "
             << Residuals<scalar>::field(mesh, "p_rgh").size() << endl;

        // The current PIMPLE pressure correction is being bypassed because the
        // initial p_rgh residual has not decreased sufficiently over
        // consecutive PIMPLE iterations.  Even though the pressure equation is
        // not solved in this iteration, append a synthetic SolverPerformance
        // entry to the p_rgh residual history.  This keeps the residual list
        // aligned with the PIMPLE iteration count, so that the skipped
        // iterations leave no gap and the final PIMPLE iteration can still be
        // reached and executed with the final relaxation factors.

        Residuals<scalar>::append(
            mesh,
            SolverPerformance<scalar>("bypass",
                                      "p_rgh",
                                      prevPimpleInitialResidual_,
                                      prevPimpleInitialResidual_,
                                      0,
                                      false,
                                      false));

        Info << "Bypass p_rgh residual history size after append = "
             << Residuals<scalar>::field(mesh, "p_rgh").size()
             << ", copied initial residual = " << prevPimpleInitialResidual_
             << endl;

        return;
    }

    if (faceMomentum)
    {
        facePressureCorrector();
    }
    else
    {
        cellPressureCorrector();
    }

    if (logPressureExtrema)
    {
        const volScalarField& p = p_;
        // const volScalarField& p = p_rgh_;  // Use this instead for p_rgh
        // extrema

        const scalarField& pCells = p.primitiveField();
        const vectorField& cellCentres = mesh.C().primitiveField();

        scalar localMinP = GREAT;
        scalar localMaxP = -GREAT;
        label localMinCell = -1;
        label localMaxCell = -1;

        forAll(pCells, celli)
        {
            if (pCells[celli] < localMinP)
            {
                localMinP = pCells[celli];
                localMinCell = celli;
            }

            if (pCells[celli] > localMaxP)
            {
                localMaxP = pCells[celli];
                localMaxCell = celli;
            }
        }

        scalar globalMinP = localMinP;
        scalar globalMaxP = localMaxP;

        reduce(globalMinP, minOp<scalar>());
        reduce(globalMaxP, maxOp<scalar>());

        const scalar minTol = small * (1 + mag(globalMinP));
        const scalar maxTol = small * (1 + mag(globalMaxP));

        label minProc =
            (localMinCell >= 0 && mag(localMinP - globalMinP) <= minTol
                 ? Pstream::myProcNo()
                 : labelMax);

        label maxProc =
            (localMaxCell >= 0 && mag(localMaxP - globalMaxP) <= maxTol
                 ? Pstream::myProcNo()
                 : labelMax);

        reduce(minProc, minOp<label>());
        reduce(maxProc, minOp<label>());

        label globalMinCell = Pstream::myProcNo() == minProc ? localMinCell : 0;
        label globalMaxCell = Pstream::myProcNo() == maxProc ? localMaxCell : 0;

        vector globalMinC = Zero;
        vector globalMaxC = Zero;

        if (Pstream::myProcNo() == minProc)
        {
            globalMinC = cellCentres[localMinCell];
        }

        if (Pstream::myProcNo() == maxProc)
        {
            globalMaxC = cellCentres[localMaxCell];
        }

        reduce(globalMinCell, sumOp<label>());
        reduce(globalMaxCell, sumOp<label>());
        reduce(globalMinC, sumOp<vector>());
        reduce(globalMaxC, sumOp<vector>());

        if (Pstream::master())
        {
            Info << "Pressure extrema for " << p.name() << ": "
                 << "min = " << globalMinP << " in cell " << globalMinCell
                 << " on processor " << minProc << " at " << globalMinC
                 << ", max = " << globalMaxP << " in cell " << globalMaxCell
                 << " on processor " << maxProc << " at " << globalMaxC << endl;
        }
    }

    fluid_.correctKinematics();
}


// ************************************************************************* //
