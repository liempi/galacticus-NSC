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
Implements a star formation rate property extractor class.
!!}

  use :: Globular_Cluster_Formation_Rates_Disks    , only : globularClusterFormationRateDisksClass
  use :: Globular_Cluster_Formation_Rates_Spheroids, only : globularClusterFormationRateSpheroidsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGlobularClusterFormationRate">
   <description>
    A node property extractor which extracts the globular cluster formation rate in a galaxy. The type of globular cluster formation rate is controlled by
    the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'' or ``{\normalfont \ttfamily total}''. The corresponding globular cluster formation
    rate is extracted as {\normalfont \ttfamily \textless\ component\textgreater\ globularClusterFormationRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorGlobularClusterFormationRate
     !!{
     A star formation rate property extractor class.
     !!}
     private
     class(globularClusterFormationRateDisksClass    ), pointer :: globularClusterFormationRateDisks_    => null()
     class(globularClusterFormationRateSpheroidsClass), pointer :: globularClusterFormationRateSpheroids_=> null()
   contains
     final     ::                globularClusterFormationRateDestructor
     procedure :: elementCount=> globularClusterFormationRateElementCount
     procedure :: extract     => globularClusterFormationRateExtract
     procedure :: names       => globularClusterFormationRateNames
     procedure :: descriptions=> globularClusterFormationRateDescriptions
     procedure :: unitsInSI   => globularClusterFormationRateUnitsInSI
  end type nodePropertyExtractorGlobularClusterFormationRate

  interface nodePropertyExtractorGlobularClusterFormationRate
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationRate} output analysis class.
     !!}
     module procedure globularClusterFormationRateConstructorParameters
     module procedure globularClusterFormationRateConstructorInternal
  end interface nodePropertyExtractorGlobularClusterFormationRate

contains

  function globularClusterFormationRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormationRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorGlobularClusterFormationRate)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters
    class(globularClusterFormationRateDisksClass           ), pointer       :: globularClusterFormationRateDisks_
    class(globularClusterFormationRateSpheroidsClass       ), pointer       :: globularClusterFormationRateSpheroids_

    !![
    <objectBuilder class="globularClusterFormationRateDisks"     name="globularClusterFormationRateDisks_"      source="parameters"/>
    <objectBuilder class="globularClusterFormationRateSpheroids" name="globularClusterFormationRateSpheroids_"  source="parameters"/>
    !!]
    self=nodePropertyExtractorGlobularClusterFormationRate(globularClusterFormationRateDisks_,globularClusterFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterFormationRateDisks_"    />
    <objectDestructor name="globularClusterFormationRateSpheroids_"/>
    !!]
    return
  end function globularClusterFormationRateConstructorParameters

  integer function globularClusterFormationRateElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterFormationRate), intent(inout) :: self
    double precision                                                   , intent(in   ) :: time
    !$GLC attributes unused :: self

    globularClusterFormationRateElementCount=2
    return
  end function globularClusterFormationRateElementCount

  function globularClusterFormationRateConstructorInternal(globularClusterFormationRateDisks_,globularClusterFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorGlobularClusterFormationRate)                        :: self
    class(globularClusterFormationRateDisksClass           ), intent(in   ), target :: globularClusterFormationRateDisks_
    class(globularClusterFormationRateSpheroidsClass       ), intent(in   ), target :: globularClusterFormationRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterFormationRateDisks_, *globularClusterFormationRateSpheroids_"/>
    !!]
    return
  end function globularClusterFormationRateConstructorInternal

  subroutine globularClusterFormationRateDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationRate} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGlobularClusterFormationRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%globularClusterFormationRateDisks_"    />
    <objectDestructor name="self%globularClusterFormationRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterFormationRateDestructor

  double precision function globularClusterFormationRateExtract(self,node,time,instance)
    !!{
    Implement a star formation rate output analysis property extractor.
    !!}
    implicit none
    double precision                                                   , dimension(:) , allocatable :: globularClusterFormationRateExtract
    class           (nodePropertyExtractorGlobularClusterFormationRate), intent(inout), target      :: self
    type            (treeNode                                         ), intent(inout), target      :: node
    double precision                                                   , intent(in   )              :: time
    type            (multiCounter                                     ), intent(inout), optional    :: instance
    !$GLC attributes unused ::self, instance
    allocate(globularClusterFormationRateExtract(2))
    globularClusterFormationRateExtract(1)=self%globularClusterFormationRateDisks_    %rate(node)
    globularClusterFormationRateExtract(2)=self%globularClusterFormationRateSpheroids_%rate(node)
    return
  end function globularClusterFormationRateExtract

  subroutine globularClusterFormationRateNames(self,time,names)
    !!{
    Return the name of the globularClusterFormationRate property.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterFormationRate), intent(inout)                             :: self
    double precision                                                   , intent(in   )                             :: time
    type            (varying_string                                   ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    allocate(names(2))
    names(1)=var_str('diskGlobularClusterFormationRate')
    names(2)=var_str('spheroidGlobularClusterFormationRate')
    return
  end subroutine globularClusterFormationRateNames

  subroutine globularClusterFormationRateDescriptions(self,time,descriptions)
    !!{
    Return a description of the globularClusterFormationRate property.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterFormationRate), intent(inout)                             :: self
    double precision                                                   , intent(in   )                             :: time
    type            (varying_string                                   ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time
    allocate(descriptions(2))
    descriptions(1)=var_str('Disk globular cluster formation rate [M☉ Gyr⁻¹].')
    descriptions(2)=var_str('Spheroidal globular cluster formation rate [M☉ Gyr⁻¹].')
    return
  end subroutine globularClusterFormationRateDescriptions

  double precision function globularClusterFormationRateUnitsInSI(self,time)
    !!{
    Return the units of the globularClusterFormationRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                                   , allocatable  , dimension(:) :: globularClusterFormationRateUnitsInSI
    class           (nodePropertyExtractorGlobularClusterFormationRate), intent(inout)               :: self
    double precision                                                   , intent(in   )               :: time
    !$GLC attributes unused :: self, time
    allocate(globularClusterFormationRateUnitsInSI(2))
    globularClusterFormationRateUnitsInSI=[massSolar/gigaYear, massSolar/gigaYear]
    return
  end function globularClusterFormationRateUnitsInSI
