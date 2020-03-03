!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !# <radiativeTransferSpectrum name="radiativeTransferSpectrumBlackBody">
  !#  <description>A photon spectrum class for blackBody spectrums.</description>
  !# </radiativeTransferSpectrum>
  type, extends(radiativeTransferSpectrumClass) :: radiativeTransferSpectrumBlackBody
     !% Implementation of a black body spectrum for radiative transfer calculations.
     private
     double precision :: temperature  , luminosityBolometric, &
          &              normalization
   contains
     procedure :: luminosity => blackBodyLuminosity
     procedure :: spectrum   => blackBodySpectrum
  end type radiativeTransferSpectrumBlackBody
  
  interface radiativeTransferSpectrumBlackBody
     !% Constructors for the {\normalfont \ttfamily blackBody} radiative transfer spectrum class.
     module procedure blackBodyConstructorParameters
     module procedure blackBodyConstructorInternal
  end interface radiativeTransferSpectrumBlackBody
  
contains

  function blackBodyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily blackBody} radiative transfer spectrum class which takes a parameter set as
    !% input.
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferSpectrumBlackBody)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: temperature, luminosityBolometric

    !# <inputParameter>
    !#   <name>temperature</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>5.0d3</defaultValue>
    !#   <description>The temperature of the black body spectrum (in Kelvin).</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>luminosityBolometric</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The bolometric luminosity of the black body spectrum (in $L_\odot$).</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=radiativeTransferSpectrumBlackBody(temperature,luminosityBolometric)
    !# <inputParametersValidate spectrum="parameters"/>
   return
  end function blackBodyConstructorParameters

  function blackBodyConstructorInternal(temperature,luminosityBolometric) result(self)
    !% Internal constructor for the {\normalfont \ttfamily blackBody} radiative transfer photon packet class.
    use :: Thermodynamics_Radiation, only : Blackbody_Radiance
    implicit none
    type            (radiativeTransferSpectrumBlackBody)                :: self
    double precision                                    , intent(in   ) :: temperature, luminosityBolometric
    !# <constructorAssign variables="temperature, luminosityBolometric"/>

    ! Compute normalization such that we get the desired bolometric luminosity when the spectrum is integrated over all
    ! wavelengths.
    self%normalization=+                   luminosityBolometric  &
         &             /Blackbody_Radiance(temperature         )
    return
  end function blackBodyConstructorInternal

  double precision function blackBodyLuminosity(self,wavelengthMinimum,wavelengthMaximum)
    !% Compute the luminosity in the given wavelength range for a black body spectrum.
    use :: FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use :: Numerical_Integration, only : Integrate    , Integrate_Done
    implicit none
    class           (radiativeTransferSpectrumBlackBody), intent(inout) :: self
    double precision                                    , intent(in   ) :: wavelengthMinimum    , wavelengthMaximum
    type            (fgsl_function                     )                :: integrandFunction
    type            (fgsl_integration_workspace        )                :: integrationWorkspace

    blackBodyLuminosity=+Integrate(                                         &
         &                                            wavelengthMinimum   , &
         &                                            wavelengthMaximum   , &
         &                                            integrand           , &
         &                                            integrandFunction   , &
         &                                            integrationWorkspace, &
         &                          toleranceAbsolute=0.0d+0              , &
         &                          toleranceRelative=1.0d-2                &
         &                         )                                        &
         &              *self%normalization
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
    
  contains

    double precision function integrand(wavelength)
      !% Integrand over black body spectrum.
      use :: Thermodynamics_Radiation, only : Blackbody_Emission, radianceTypeWavelength
      implicit none
      double precision, intent(in   ) :: wavelength

      integrand=Blackbody_Emission(wavelength,self%temperature,radianceTypeWavelength)
      return
    end function integrand

  end function blackBodyLuminosity

  double precision function blackBodySpectrum(self,wavelength)
    !% Return the spectrum of the black body.
    use :: Thermodynamics_Radiation, only : Blackbody_Emission, radianceTypeWavelength
    implicit none
    class           (radiativeTransferSpectrumBlackBody), intent(inout) :: self
    double precision                                    , intent(in   ) :: wavelength

    blackBodySpectrum=+                              self%normalization                         &
         &            *Blackbody_Emission(wavelength,self%temperature  ,radianceTypeWavelength)
    return
  end function blackBodySpectrum