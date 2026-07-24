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

#include "GidaspowConductivity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace conductivityModels
{
defineTypeNameAndDebug(Gidaspow, 0);

addToRunTimeSelectionTable(conductivityModel, Gidaspow, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Gidaspow::Gidaspow(
    const dictionary& coeffDict)
: conductivityModel(coeffDict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::conductivityModels::Gidaspow::~Gidaspow() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::conductivityModels::Gidaspow::kappa(
    const volScalarField& alpha1,
    const volScalarField& Theta,
    const volScalarField& g0,
    const volScalarField& sumAlphaGs0,
    const volScalarField& beta,
    const volScalarField& rho1,
    const volScalarField& da,
    const dimensionedScalar& e) const
{
    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const scalar Pi = constant::mathematical::pi;

    // beta is kept in the function signature for consistency with the common
    // conductivity-model interface, but it is intentionally not used here.
    (void)beta;

    // g0 is kept in the function signature for consistency with the common
    // conductivity-model interface. In the present split-consistent
    // polydisperse form it is not used explicitly, because the monodisperse
    // factor 1/g0 is recovered from alpha1/sumAlphaGs0 when
    //
    //     sumAlphaGs0 -> alpha1*g0.
    (void)g0;

    // eta = (1 + e)/2
    // This is the standard compact notation used to rewrite the original
    // monodisperse Gidaspow conductivity in modular form.
    const dimensionedScalar eta = 0.5 * (1.0 + e);

    // Base Gidaspow transport prefactor in modular form:
    //
    //     kappaGid** = 75*rho*d*sqrt(pi*Theta)/(384*eta)
    //
    // In the monodisperse Gidaspow expression, the transport prefactor appears
    // as kappaGid**/g0 and the crowding factor is alpha*g0:
    //
    //     kappa =
    //         (kappaGid**/g0)
    //         * [ (1 + 12/5*eta*alpha*g0)^2
    //             + 512/(25*pi)*eta^2*(alpha*g0)^2 ].
    //
    // For the multi-solid formulation we use
    //
    //     alpha*g0 -> sum_j alpha_j*g0_ij = sumAlphaGs0.
    //
    // However, unlike the shear viscosity in divDevTau, the granular
    // conductivity kappa enters the granular-temperature equation directly as
    //
    //     div(kappa grad(Theta)),
    //
    // without an additional external factor alpha_i. Therefore kappa itself
    // must carry the abundance of the current solid phase.
    //
    // To preserve both the monodisperse limit and split-phase consistency, the
    // monodisperse factor 1/g0 is rewritten as
    //
    //     alpha_i/(alpha_i*g0)
    //
    // and generalized as
    //
    //     alpha_i/sumAlphaGs0.
    //
    // This gives the original monodisperse expression when
    // sumAlphaGs0 = alpha_i*g0, but, if one solid phase is split into two
    // identical phases, each split phase receives a conductivity proportional
    // to its own alpha_i. The sum of the two split conductivities then
    // reconstructs the unsplit conductivity.
    const volScalarField kappaStarStar((75.0 * rho1 * da * sqrtPi * sqrt(Theta))
                                       / (384.0 * eta));

    // Split-consistent polydisperse Gidaspow conductivity without drag
    // correction:
    //
    //     kappa =
    //         kappaGid** * alpha_i/sumAlphaGs0
    //         * [ (1 + 12/5*eta*sumAlphaGs0)^2
    //             + 512/(25*pi)*eta^2*(sumAlphaGs0)^2 ].
    //
    // In expanded form, the terms proportional to (sumAlphaGs0)^2 are therefore
    // effectively proportional to alpha_i*sumAlphaGs0, not to
    // (sumAlphaGs0)^2/g0. This is the required scaling for split consistency.
    return kappaStarStar * alpha1 / sumAlphaGs0
         * (sqr(1.0 + 12.0 / 5.0 * eta * sumAlphaGs0)
            + 512.0 / (25.0 * Pi) * sqr(eta) * sqr(sumAlphaGs0));
}

// ************************************************************************* //
