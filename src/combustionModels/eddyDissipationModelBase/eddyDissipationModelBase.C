/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd
-------------------------------------------------------------------------------
License

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

#include "eddyDissipationModelBase.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationModelBase<ReactionThermo, ThermoType>::eddyDissipationModelBase
(
    const word& modelType,
    const ReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion<ReactionThermo, ThermoType>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    CEDC_(this->coeffs().template lookup <scalar>("CEDC")),
    finiteRateModel_(this->coeffs().template lookupOrDefault <Switch>("FiniteRateModel",false)),
    A_(0.0),
    beta_(0.0),
    Ta_(0.0),
    n1_(0.0),
    n2_(0.0),
    wFR_
    (
        IOobject
        (
            "wFR",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    ),
    wEDM_
    (
        IOobject
        (
            "wEDM",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimMass/dimVolume/dimTime, 0)
    )
{
    // By default the FR/EDM model is off
    if (finiteRateModel_)
    {
        A_ = this->coeffs().template lookup <scalar>("A");
        beta_ = this->coeffs().template lookup <scalar>("beta");
        Ta_ = this->coeffs().template lookup <scalar>("Ta");
        n1_ = this->coeffs().template lookup <scalar>("n1");
        n2_ = this->coeffs().template lookup <scalar>("n2");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationModelBase<ReactionThermo, ThermoType>::
~eddyDissipationModelBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationModelBase<ReactionThermo, ThermoType>::rtTurb() const
{
    return
        CEDC_*this->turbulence().epsilon()
      / max
        (
            this->turbulence().k(),
            dimensionedScalar("SMALL",sqr(dimVelocity), SMALL)
        );
}


template<class ReactionThermo, class ThermoType>
void eddyDissipationModelBase<ReactionThermo, ThermoType>::correct()
{
    this->wFuel_ == dimensionedScalar(dimMass/dimVolume/dimTime, Zero);

    if (true)
    {
        this->fresCorrect();

        const label fuelI = this->fuelIndex();

        const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];

        const dimensionedScalar s = this->s();

        if (this->thermo_.composition().contains("O2"))
        {
            const volScalarField& YO2 = this->thermo().composition().Y("O2");

            // Eddy dissipation model reaction rate
            wEDM_ == this->rho() * min(YFuel, YO2/s.value()) * timeScale();
            
            if (finiteRateModel_)
            {
                // Dimension transformation
                dimensionedScalar _Ta(dimTemperature,Ta_);
                dimensionedScalar _A(dimVolume/dimMass/dimTime,A_);
                // Arrhenius reaction rate
                wFR_ = sqr(this->rho())*_A*pow(YFuel,n1_)*pow(YO2,n2_)*exp(-_Ta/this->thermo().T());
                // Combined Finite Rate /Eddy Dissipation model
                this->wFuel_ == min(wFR_, wEDM_);
                Info<<"wFR min/max="<<min(wFR_).value()<<","<<max(wFR_).value()<<endl;
                Info<<"wEDM min/max="<<min(wEDM_).value()<<","<<max(wEDM_).value()<<endl;
            }
            else
            {
                // Eddy Dissipation model
                this->wFuel_ == wEDM_;
            }

        }
        else
        {
            FatalErrorInFunction
                << "You selected a combustion model that requires O2 mass"
                << " to be present in the mixture"
                << exit(FatalError);
        }
    }
}


template<class ReactionThermo, class ThermoType>
bool eddyDissipationModelBase<ReactionThermo, ThermoType>::read()
{
    if (singleStepCombustion<ReactionThermo, ThermoType>::read())
    {
        this->coeffs().lookup("CEDC") >> CEDC_;
        this->coeffs().lookup("FiniteRateModel") >> finiteRateModel_;
        if (finiteRateModel_)
        {
            this->coeffs().lookup("A") >> A_;
            this->coeffs().lookup("beta") >> beta_;
            this->coeffs().lookup("Ta") >> Ta_;
            this->coeffs().lookup("n1") >> n1_;
            this->coeffs().lookup("n2") >> n2_;
        }
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
