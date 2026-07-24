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

#include "LunPressure.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
defineTypeNameAndDebug(Lun, 0);

addToRunTimeSelectionTable(granularPressureModel, Lun, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Lun::Lun(
    const dictionary& coeffDict)
: granularPressureModel(coeffDict)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Lun::~Lun() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::granularPressureCoeff(
    const phaseModel& phase1,
    const phaseModel& continuousPhase,
    const PtrList<volScalarField>& g0_im,
    const volScalarField& rho1,
    const dimensionedScalar& e) const
{
    // pressure coeff from Eq. B3 MFIX Equations 2012

    volScalarField alpha1 = phase1;
    const phaseSystem& fluid = phase1.fluid();

    volScalarField pCoeff(alpha1 * rho1);
    dimensionedScalar eta = 0.5 * (1.0 + e);

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];

        if (&phase != &continuousPhase)
        {
            const volScalarField& alpha = phase;

            pCoeff += alpha1 * rho1 * (4.0 * eta * alpha * g0_im[phasei]);
        }
    }

    return 1.0 * pCoeff;
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::
    granularPressureCoeffPrime(const phaseModel& phase1,
                               const phaseModel& continuousPhase,
                               const PtrList<volScalarField>& g0_im,
                               const PtrList<volScalarField>& g0prime_im,
                               const volScalarField& rho1,
                               const dimensionedScalar& e) const
{
    const phaseSystem& fluid = phase1.fluid();

    // eta = (1 + e)/2, used to write the collisional contribution
    // in the compact form 4*eta = 2*(1 + e).
    const dimensionedScalar eta("eta", dimless, 0.5 * (1.0 + e.value()));

    // This function returns the derivative of the granular-pressure
    // coefficient with respect to alpha_1, i.e. the coefficient C'_1
    // such that:
    //
    //     p'_1 = Theta_1 * C'_1 + p'_{1,fric}
    //
    // For the Lun model:
    //
    //     p_{1,k} = rho_1 * alpha_1 * [1 + 2(1 + e) * Sigma_1] * Theta_1
    //
    // with
    //
    //     Sigma_1 = sum_m alpha_m * g0_{1m}.
    //
    // Therefore, in the frozen-g0 approximation:
    //
    //     C'_1 = rho_1 + 2(1 + e) * rho_1 * Sigma_1
    //          = rho_1 + 4*eta*rho_1*Sigma_1.
    //
    // The first term below is the kinetic contribution "rho_1".
    tmp<volScalarField> tpCoeffPrime(
        volScalarField::New("granularPressureCoeffPrime", rho1));

    volScalarField& pCoeffPrime = tpCoeffPrime.ref();

    // Add the collisional contribution in frozen-g0 form:
    //
    //     4*eta*rho_1*Sigma_1
    //   = 4*eta*rho_1*sum_m(alpha_m*g0_{1m})
    //
    // This is the split-consistent multi-solid form used by default.
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase2 = fluid.phases()[phasei];

        if (&phase2 != &continuousPhase)
        {
            const volScalarField& alpha2 = phase2;

            pCoeffPrime += rho1 * (4.0 * eta * alpha2 * g0_im[phasei]);
        }
    }

    // Optional correction:
    // when includeG0primeInPPrime() is enabled, add the derivative terms
    // associated with the alpha_1-dependence of alpha_1*g0_{11}.
    //
    // In the monodisperse limit:
    //
    //     d/dalpha_1 [alpha_1^2 g0_{11}]
    //   = 2*alpha_1*g0_{11} + alpha_1^2*g0'_{11}.
    //
    // One factor alpha_1*g0_{11} is already contained in the frozen-g0 sum
    // above. The additional term added here is therefore:
    //
    //     alpha_1*g0_{11} + alpha_1^2*g0'_{11},
    //
    // so that the total result recovers the full analytical derivative of
    // the original monodisperse Lun expression.
    if (fluid.includeG0primeInPPrime())
    {
        const volScalarField& alpha1 = phase1;
        const label index1 = phase1.index();

        pCoeffPrime +=
            rho1
            * (4.0 * eta
               * (alpha1 * g0_im[index1] + sqr(alpha1) * g0prime_im[index1]));
    }

    return tpCoeffPrime;
}

// ************************************************************************* //
