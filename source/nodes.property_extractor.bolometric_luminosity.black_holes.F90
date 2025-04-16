!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!!{
Contains a module which implements an AGN bolometric luminosity property extractor class.
!!}
  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorAGNBolometricLuminosity">
   <description>An AGN bolometric luminosity property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorAGNBolometricLuminosity
     !!{
     An AGNBolometricLuminosity property extractor class.
     !!}
     private
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class(accretionDisksClass        ), pointer :: accretionDisks_         => null()
   contains
     final     ::                agnBolometricLuminosityDestructor
     procedure :: extract     => agnBolometricLuminosityExtract
     procedure :: name        => agnBolometricLuminosityName
     procedure :: description => agnBolometricLuminosityDescription
     procedure :: unitsInSI   => agnBolometricLuminosityUnitsInSI
  end type nodePropertyExtractorAGNBolometricLuminosity

  interface nodePropertyExtractorAGNBolometricLuminosity
     !!{
     Constructors for the ``AGNBolometricLuminosity'' output analysis class.
     !!}
     module procedure agnBolometricLuminosityConstructorParameters
     module procedure agnBolometricLuminosityConstructorInternal
  end interface nodePropertyExtractorAGNBolometricLuminosity

contains

  function agnBolometricLuminosityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily AGNBolometricLuminosity} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorAGNBolometricLuminosity)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(blackHoleAccretionRateClass                 ), pointer       :: blackHoleAccretionRate_
    class(accretionDisksClass                         ), pointer       :: accretionDisks_

    !![
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_"  source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"          source="parameters"/> 
    !!]
    self=nodePropertyExtractorAGNBolometricLuminosity(blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_" />
    <objectDestructor name="accretionDisks_"/>
    !!]
    return
  end function agnBolometricLuminosityConstructorParameters

  function agnBolometricLuminosityConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily AGNBolometricLuminosity} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorAGNBolometricLuminosity)                        :: self
    class(blackHoleAccretionRateClass                 ), intent(in   ), target :: blackHoleAccretionRate_
    class(accretionDisksClass                         ), intent(in   ), target :: accretionDisks_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    return
  end function agnBolometricLuminosityConstructorInternal

  subroutine agnBolometricLuminosityDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily AGNBolometricLuminosity} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorAGNBolometricLuminosity), intent(inout) :: self
    !![
    <objectDestructor name="self%blackHoleAccretionRate_" />
    <objectDestructor name="self%accretionDisks_"/>
    !!]
    return
  end subroutine agnBolometricLuminosityDestructor

  double precision function agnBolometricLuminosityExtract(self,node,instance)
    !!{
    Implement an AGN bolometric luminosity extractor.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gigaYear              , luminositySolar, massSolar
    use :: Numerical_Constants_Physical    , only : speedLight
    implicit none
    class           (nodePropertyExtractorAGNBolometricLuminosity), intent(inout), target   :: self
    type            (treeNode                                    ), intent(inout), target   :: node
    type            (multiCounter                                ), intent(inout), optional :: instance
    class           (nodeComponentBlackHole                      )               , pointer  :: blackHole
    double precision                                                                        :: rateAccretionSpheroid, rateAccretionHotHalo           , &
         &                                                                                     rateAccretion        , rateAccretionNuclearStarCluster, &
         &                                                                                     efficiencyRadiative
    !$GLC attributes unused :: instance
    
    blackHole                      => node%blackHole()
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion                  = +                    rateAccretionSpheroid                                                         &
         &                           +                    rateAccretionHotHalo                                                          &
         &                           +                    rateAccretionNuclearStarCluster
    efficiencyRadiative            = self%accretionDisks_%efficiencyRadiative  (blackHole,rateAccretion                               )
    agnBolometricLuminosityExtract = +efficiencyRadiative    &
          &                          *rateAccretion          &
          &                          *massSolar              &
          &                          *speedLight         **2 &
          &                          /gigaYear               &
          &                          /luminositySolar
    return
  end function agnBolometricLuminosityExtract

  function agnBolometricLuminosityName(self)
    !!{
    Return the names of the {\normalfont \ttfamily AGNBolometricLuminosity} properties.
    !!}
    implicit none
    type (varying_string                              )                :: agnBolometricLuminosityName
    class(nodePropertyExtractorAGNBolometricLuminosity), intent(inout) :: self
    !$GLC attributes unused :: self
    
    agnBolometricLuminosityName=var_str('AGNBolometricLuminosity')
    return
  end function agnBolometricLuminosityName

  function agnBolometricLuminosityDescription(self)
    !!{
    Return descriptions of the {\normalfont \ttfamily AGNBolometricLuminosity} properties.
    !!}
    implicit none
    type (varying_string                              )                :: agnBolometricLuminosityDescription
    class(nodePropertyExtractorAGNBolometricLuminosity), intent(inout) :: self
    !$GLC attributes unused :: self
    
    agnBolometricLuminosityDescription=var_str('AGN bolometric luminosity [ergs/s]')
    return
  end function agnBolometricLuminosityDescription

  double precision function agnBolometricLuminosityUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily AGNBolometricLuminosity} properties in the SI system.
    !!}
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    class           (nodePropertyExtractorAGNBolometricLuminosity), intent(inout) :: self
    !$GLC attributes unused :: self
    
    agnBolometricLuminosityUnitsInSI=ergs
    return
  end function agnBolometricLuminosityUnitsInSI

