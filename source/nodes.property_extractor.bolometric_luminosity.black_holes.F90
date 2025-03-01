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
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorAGNBolometricLuminosity
     !!{
     An AGNBolometricLuminosity property extractor class.
     !!}
     private
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class(accretionDisksClass        ), pointer :: accretionDisks_         => null()
   contains
     final     ::                 AGNBolometricLuminosityDestructor
     procedure :: elementCount => AGNBolometricLuminosityElementCount
     procedure :: extract      => AGNBolometricLuminosityExtract
     procedure :: names        => AGNBolometricLuminosityNames
     procedure :: descriptions => AGNBolometricLuminosityDescriptions
     procedure :: unitsInSI    => AGNBolometricLuminosityUnitsInSI
  end type nodePropertyExtractorAGNBolometricLuminosity

  interface nodePropertyExtractorAGNBolometricLuminosity
     !!{
     Constructors for the ``AGNBolometricLuminosity'' output analysis class.
     !!}
     module procedure AGNBolometricLuminosityConstructorParameters
     module procedure AGNBolometricLuminosityConstructorInternal
  end interface nodePropertyExtractorAGNBolometricLuminosity

contains

  function AGNBolometricLuminosityConstructorParameters(parameters) result(self)
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
  end function AGNBolometricLuminosityConstructorParameters

  function AGNBolometricLuminosityConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
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
  end function AGNBolometricLuminosityConstructorInternal

  subroutine AGNBolometricLuminosityDestructor(self)
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
  end subroutine AGNBolometricLuminosityDestructor

  integer function AGNBolometricLuminosityElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorAGNBolometricLuminosity), intent(inout) :: self

    AGNBolometricLuminosityElementCount=1
    return
  end function AGNBolometricLuminosityElementCount

  function AGNBolometricLuminosityExtract(self,node,instance) result(AGNBolometricLuminosity)
    !!{
    Implement an ICM X-ray properties extractor.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gigaYear              , luminositySolar, massSolar
    use :: Numerical_Constants_Physical    , only : speedLight
    implicit none
    double precision                                              , dimension(:,:), allocatable :: AGNBolometricLuminosity
    class           (nodePropertyExtractorAGNBolometricLuminosity), intent(inout)               :: self
    type            (treeNode                                    ), intent(inout)               :: node
    type            (multiCounter                                ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole                      )                , pointer     :: blackHole
    double precision                                                                            :: rateAccretionSpheroid  , rateAccretionHotHalo           , &
         &                                                                                         rateAccretion          , rateAccretionNuclearStarCluster, &
         &                                                                                         efficiencyRadiative
    integer                                                                                     :: i                      , countBlackHoles
    !$GLC attributes unused :: instance
    
    countBlackHoles=node%blackHoleCount()
    allocate(AGNBolometricLuminosity(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole                   => node%blackHole(instance=i)
       call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
       rateAccretion               = +                    rateAccretionSpheroid                                                         &
            &                        +                    rateAccretionHotHalo                                                          &
            &                        +                    rateAccretionNuclearStarCluster
       efficiencyRadiative         = self%accretionDisks_%efficiencyRadiative  (blackHole,rateAccretion                               )
       AGNBolometricLuminosity(i,1)= +efficiencyRadiative    &
            &                        *rateAccretion          &
            &                        *massSolar              &
            &                        *speedLight         **2 &
            &                        /gigaYear               &
            &                        /luminositySolar
    end do
    return
  end function AGNBolometricLuminosityExtract

  subroutine AGNBolometricLuminosityNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily AGNBolometricLuminosity} properties.
    !!}
    implicit none
    class(nodePropertyExtractorAGNBolometricLuminosity), intent(inout)                             :: self
    type (varying_string                              ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self
    
    allocate(names(1))
    names(1)=var_str('AGNBolometricLuminosity')
    return
  end subroutine AGNBolometricLuminosityNames

  subroutine AGNBolometricLuminosityDescriptions(self,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily AGNBolometricLuminosity} properties.
    !!}
    implicit none
    class(nodePropertyExtractorAGNBolometricLuminosity), intent(inout)                             :: self
    type (varying_string                              ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self
    
    allocate(descriptions(1))
    descriptions(1)=var_str('AGN luminosity of the ICM [ergs/s]')
    return
  end subroutine AGNBolometricLuminosityDescriptions

  function AGNBolometricLuminosityUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily AGNBolometricLuminosity} properties in the SI system.
    !!}
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    double precision                                              , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorAGNBolometricLuminosity), intent(inout)              :: self
    !$GLC attributes unused :: self
    
    allocate(unitsInSI(1))
    unitsInSI(1)=ergs
    return
  end function AGNBolometricLuminosityUnitsInSI

