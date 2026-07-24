/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenPDAC.
    This file was derived from the multiphaseEuler solver in OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "CarnahanStarling.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radialModels
{
defineTypeNameAndDebug(CarnahanStarling, 0);

addToRunTimeSelectionTable(radialModel, CarnahanStarling, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radialModels::CarnahanStarling::CarnahanStarling(const dictionary& dict,
                                                       const phaseSystem& fluid)
: radialModel(dict, fluid)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radialModels::CarnahanStarling::~CarnahanStarling() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField> Foam::radialModels::CarnahanStarling::g0(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax) const
{
    const volScalarField& alphai = phasei;
    const label& indexi = phasei.index();
    const phaseSystem& fluid = phasei.fluid();

    PtrList<volScalarField> g0_im(fluid.phases().size());

    volScalarField alphas = alphai;
    volScalarField eta2(alphai / phasei.d());

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];
        if ((&phase != &continuousPhase) and !(phaseIdx == indexi))
        {
            alphas += phase;
            eta2 += phase / phase.d();
        }
    }

    const volScalarField denominatorTerm(1.0 - alphas);

    forAll(g0_im, iter)
    {
        const phaseModel& phasem = fluid.phases()[iter];

        if (&phasem != &continuousPhase)
        {
            const volScalarField di(phasei.d());
            const volScalarField dm(phasem.d());
            const dimensionedScalar smallD("smallD", dimLength, ROOTVSMALL);
            volScalarField term_d(di * dm / (di + dm + smallD));

            g0_im.set(iter,
                      new volScalarField(
                          "g0_im" + phasei.name() + "_" + phasem.name(),
                          1.0 / denominatorTerm
                              + 3.0 * term_d * eta2 / sqr(denominatorTerm)
                              + 2.0 * sqr(term_d) * sqr(eta2)
                                    / pow3(denominatorTerm)));
        }
    }

    return g0_im;
}

Foam::PtrList<Foam::volScalarField>
Foam::radialModels::CarnahanStarling::g0prime(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax) const
{
    // Note:
    // In the current OpenPDAC design, this function returns
    //
    //     g0prime_im[m] = d g0_{i,m} / d alpha_i
    //
    // i.e. the derivative of all pair radial-distribution functions involving
    // phase i with respect to the volume fraction of the same phase i.
    //
    // For the Carnahan-Starling model implemented here:
    //
    //   g0_{i,m} =
    //       1/D
    //     + 3*A_{i,m}*eta2/D^2
    //     + 2*A_{i,m}^2*eta2^2/D^3
    //
    // where
    //
    //   D    = 1 - alpha_s,tot
    //   eta2 = sum_j(alpha_j/d_j)
    //   A_{i,m} = d_i d_m / (d_i + d_m)
    //
    // Since alpha_s,tot = sum_j alpha_j over solids, one has
    //
    //   dD/d(alpha_i) = -1
    //   d(1/D)/d(alpha_i) = 1/D^2
    //
    // and, crucially,
    //
    //   d(eta2)/d(alpha_i) = 1/d_i
    //
    // not eta2/alpha_s,tot.
    //
    // The previous implementation used eta2/alphas, which is only correct in
    // special cases and is generally wrong for polydisperse mixtures.

    const phaseSystem& fluid = phasei.fluid();
    const label indexi = phasei.index();

    PtrList<volScalarField> g0prime_im(fluid.phases().size());

    // Total solids volume fraction: alpha_s,tot = sum_j alpha_j
    volScalarField alphas = phasei;

    // Mixture moment: eta2 = sum_j(alpha_j/d_j)
    volScalarField eta2(phasei / phasei.d());

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if ((&phase != &continuousPhase) && (phaseIdx != indexi))
        {
            alphas += phase;
            eta2 += phase / phase.d();
        }
    }

    // Denominator appearing in the Carnahan-Starling expression
    const volScalarField denominatorTerm(1.0 - alphas);

    // Exact derivative of eta2 with respect to alpha_i:
    //
    //   d/d(alpha_i) [ sum_j(alpha_j/d_j) ] = 1/d_i
    //
    const volScalarField dEta2dAlphai(1.0 / phasei.d());

    forAll(g0prime_im, iter)
    {
        const phaseModel& phasem = fluid.phases()[iter];

        if (&phasem != &continuousPhase)
        {
            const volScalarField di(phasei.d());
            const volScalarField dm(phasem.d());

            // Pair size factor:
            //
            //   A_{i,m} = d_i d_m / (d_i + d_m)
            //
            // A small denominator safeguard is kept for robustness.
            const dimensionedScalar smallD("smallD", dimLength, ROOTVSMALL);
            const volScalarField Aim(di * dm / (di + dm + smallD));

            // Derivative of:
            //
            //   T1 = 1/D
            //   T2 = 3*A*eta2/D^2
            //   T3 = 2*A^2*eta2^2/D^3
            //
            // with respect to alpha_i.
            //
            // Since dD/d(alpha_i) = -1:
            //
            //   d(T1) = 1/D^2
            //   d(T2) = 3*A * [ dEta2/D^2 + 2*eta2/D^3 ]
            //   d(T3) = 2*A^2 * [ 2*eta2*dEta2/D^3 + 3*eta2^2/D^4 ]
            //
            tmp<volScalarField> dT1 = 1.0 / sqr(denominatorTerm);

            tmp<volScalarField> dT2 = 3.0 * Aim
                                    * (dEta2dAlphai / sqr(denominatorTerm)
                                       + 2.0 * eta2 / pow3(denominatorTerm));

            tmp<volScalarField> dT3 =
                2.0 * sqr(Aim)
                * (2.0 * eta2 * dEta2dAlphai / pow3(denominatorTerm)
                   + 3.0 * sqr(eta2) / pow(denominatorTerm, 4.0));

            g0prime_im.set(iter,
                           new volScalarField("g0prime_im" + phasei.name() + "_"
                                                  + phasem.name(),
                                              dT1 + dT2 + dT3));
        }
    }

    return g0prime_im;
}

// ************************************************************************* //
