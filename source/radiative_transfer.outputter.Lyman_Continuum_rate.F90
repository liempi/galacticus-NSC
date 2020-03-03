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

  !# <radiativeTransferOutputter name="radiativeTransferOutputterLymanContinuumRate">
  !#  <description>A radiative transfer outputter class which outputs the Lyman continuum photon emission rate.</description>
  !# </radiativeTransferOutputter>
  type, extends(radiativeTransferOutputterClass) :: radiativeTransferOutputterLymanContinuumRate
     !% Implementation of a radiative transfer outputter class which outputs the Lyman continuum photon emission rate.
     private
     double precision :: lymanContinuumRateEscaping
   contains
     procedure :: reset               => lymanContinuumRateReset
     procedure :: sourceProperties    => lymanContinuumRateSourceProperties
     procedure :: photonPacketEscapes => lymanContinuumRatePhotonPacketEscapes
     procedure :: finalize            => lymanContinuumRateFinalize
     procedure :: output              => lymanContinuumRateOutput
  end type radiativeTransferOutputterLymanContinuumRate

  interface radiativeTransferOutputterLymanContinuumRate
     !% Constructors for the {\normalfont \ttfamily lymanContinuumRate} radiative transfer outputter packet class.
     module procedure lymanContinuumRateConstructorParameters
     module procedure lymanContinuumRateConstructorInternal
  end interface radiativeTransferOutputterLymanContinuumRate
  
contains

  function lymanContinuumRateConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily lymanContinuumRate} radiative transfer outputter class which takes a parameter set as
    !% input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(radiativeTransferOutputterLymanContinuumRate)                :: self
    type(inputParameters                             ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=radiativeTransferOutputterLymanContinuumRate()
    return
  end function lymanContinuumRateConstructorParameters
  
  function lymanContinuumRateConstructorInternal() result(self)
    !% Internal constructor for the {\normalfont \ttfamily lymanContinuumRate} radiative transfer outputter class.
    implicit none
    type(radiativeTransferOutputterLymanContinuumRate) :: self
    
    self%lymanContinuumRateEscaping=0.0d0
    return
  end function lymanContinuumRateConstructorInternal
  
  subroutine lymanContinuumRateReset(self)
    !% Reset the accumulated Lyman continuum photon escape rate.
    implicit none
    class(radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self

    self%lymanContinuumRateEscaping=0.0d0
    return
  end subroutine lymanContinuumRateReset

  subroutine lymanContinuumRateSourceProperties(self,radiativeTransferSource_,outputGroup)
    !% Compute and output the Lyman continuum photon emission rate.
    use :: FGSL                      , only : fgsl_function                     , fgsl_integration_workspace
    use :: IO_HDF5                   , only : hdf5Access
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen
    use :: Numerical_Integration     , only : Integrate                         , Integrate_Done
    implicit none
    class           (radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self
    class           (radiativeTransferSourceClass                ), intent(inout) :: radiativeTransferSource_
    type            (hdf5Object                                  ), intent(inout) :: outputGroup
    type            (fgsl_function                               )                :: integrandFunction
    type            (fgsl_integration_workspace                  )                :: integrationWorkspace
    double precision                                                              :: rateLymanContinuum
    !GCC$ attributes unused :: self
    
    rateLymanContinuum=+Integrate(                                         &
         &                                            1.0d-6*lymanSeriesLimitWavelengthHydrogen   , &
         &                                            lymanSeriesLimitWavelengthHydrogen   , &
         &                                            integrand           , &
         &                                            integrandFunction   , &
         &                                            integrationWorkspace, &
         &                          toleranceAbsolute=0.0d+0              , &
         &                          toleranceRelative=1.0d-2                &
         &                         )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    !$ call hdf5Access%set  ()
    call outputGroup%writeAttribute(rateLymanContinuum,'rateLymanContinuumEmitted')
    !$ call hdf5Access%unset()
    return
    
  contains

    double precision function integrand(wavelength)
      !% Integrand over the source spectrum.
      use :: Numerical_Constants_Physical    , only : plancksConstant  , speedLight
      use :: Numerical_Constants_Units       , only : angstromsPerMeter
      use :: Numerical_Constants_Astronomical, only : luminositySolar
      implicit none
      double precision, intent(in   ) :: wavelength
      double precision                :: energyPhoton
      
      energyPhoton=+plancksConstant                               &
           &       *speedLight                                    &
           &       *angstromsPerMeter                             &
           &       /wavelength
      integrand   =+radiativeTransferSource_%spectrum(wavelength) &
           &       *luminositySolar                               &
           &       /energyPhoton
      return
    end function integrand

  end subroutine lymanContinuumRateSourceProperties

  subroutine lymanContinuumRatePhotonPacketEscapes(self,photonPacket)
    !% Process an escaping photon packet.
    use :: Numerical_Constants_Atomic      , only : lymanSeriesLimitWavelengthHydrogen
    use :: Numerical_Constants_Physical    , only : plancksConstant                   , speedLight
    use :: Numerical_Constants_Units       , only : angstromsPerMeter
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    class           (radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass          ), intent(inout) :: photonPacket
    double precision                                                              :: energyPhoton

    if (photonPacket%wavelength() < lymanSeriesLimitWavelengthHydrogen) then
       energyPhoton                   =+plancksConstant                           &
            &                          *speedLight                                &
            &                          *angstromsPerMeter                         &
            &                          /photonPacket%wavelength                ()
       self%lymanContinuumRateEscaping=+self        %lymanContinuumRateEscaping   &
            &                          +photonPacket%luminosity                () &
            &                          *luminositySolar                           &
            &                          /energyPhoton
    end if
    return
  end subroutine lymanContinuumRatePhotonPacketEscapes

  subroutine lymanContinuumRateFinalize(self)
    !% Finalize the Lyman continuum photon escape rate.
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class(radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self

    ! Sum the Lyc rate across all MPI processes.
    self%lymanContinuumRateEscaping=mpiSelf%sum(self%lymanContinuumRateEscaping)
    return
  end subroutine lymanContinuumRateFinalize

  subroutine lymanContinuumRateOutput(self,outputGroup)
    !% Output the Lyman continuum photon escape rate.
    use :: IO_HDF5, only : hdf5Access
    implicit none
    class(radiativeTransferOutputterLymanContinuumRate), intent(inout) :: self
    type (hdf5Object                                  ), intent(inout) :: outputGroup

    !$ call hdf5Access%set  ()
    call outputGroup%writeAttribute(self%lymanContinuumRateEscaping,'rateLymanContinuumEscaping')
    !$ call hdf5Access%unset()
    return
  end subroutine lymanContinuumRateOutput