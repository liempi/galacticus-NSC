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
Implements a star Dissolution rate property extractor class.
!!}

  use :: Globular_Cluster_Dissolution_Rates_Disks    , only : globularClusterDissolutionRateDisksClass
  use :: Globular_Cluster_Dissolution_Rates_Spheroids, only : globularClusterDissolutionRateSpheroidsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGlobularClusterDissolutionRate">
   <description>
    A node property extractor which extracts the globular cluster dissolution rate in a galaxy. The type of globular dissolution rate is controlled by
    the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'' or ``{\normalfont \ttfamily total}''. The corresponding globular cluster dissolution
    rate is extracted as {\normalfont \ttfamily \textless\ component\textgreater\ globularClusterDissolutionRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorGlobularClusterDissolutionRate
     !!{
     A star Dissolution rate property extractor class.
     !!}
     private
     class(globularClusterDissolutionRateDisksClass    ), pointer :: globularClusterDissolutionRateDisks_    => null()
     class(globularClusterDissolutionRateSpheroidsClass), pointer :: globularClusterDissolutionRateSpheroids_=> null()
   contains
     final     ::                 globularClusterDissolutionRateDestructor
     procedure :: elementCount => globularClusterDissolutionRateElementCount
     procedure :: extract      => globularClusterDissolutionRateExtract
     procedure :: names        => globularClusterDissolutionRateNames
     procedure :: descriptions => globularClusterDissolutionRateDescriptions
     procedure :: unitsInSI    => globularClusterDissolutionRateUnitsInSI
  end type nodePropertyExtractorGlobularClusterDissolutionRate

  interface nodePropertyExtractorGlobularClusterDissolutionRate
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterDissolutionRate} output analysis class.
     !!}
     module procedure globularClusterDissolutionRateConstructorParameters
     module procedure globularClusterDissolutionRateConstructorInternal
  end interface nodePropertyExtractorGlobularClusterDissolutionRate

contains

  function globularClusterDissolutionRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterDissolutionRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorGlobularClusterDissolutionRate)                :: self
    type (inputParameters                                    ), intent(inout) :: parameters
    class(globularClusterDissolutionRateDisksClass           ), pointer       :: globularClusterDissolutionRateDisks_
    class(globularClusterDissolutionRateSpheroidsClass       ), pointer       :: globularClusterDissolutionRateSpheroids_

    !![
    <objectBuilder class="globularClusterDissolutionRateDisks"     name="globularClusterDissolutionRateDisks_"     source="parameters"/>
    <objectBuilder class="globularClusterDissolutionRateSpheroids" name="globularClusterDissolutionRateSpheroids_" source="parameters"/>
    !!]
 
    self=nodePropertyExtractorGlobularClusterDissolutionRate(globularClusterDissolutionRateDisks_,globularClusterDissolutionRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterDissolutionRateDisks_"    />
    <objectDestructor name="globularClusterDissolutionRateSpheroids_"/>
    !!]
    return
  end function globularClusterDissolutionRateConstructorParameters

  integer function globularClusterDissolutionRateElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout) :: self
    double precision                                                     , intent(in   ) :: time
    !$GLC attributes unused :: self

    globularClusterDissolutionRateElementCount=2
    return
  end function globularClusterDissolutionRateElementCount

  function globularClusterDissolutionRateConstructorInternal(globularClusterDissolutionRateDisks_,globularClusterDissolutionRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterDissolutionRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorGlobularClusterDissolutionRate)                        :: self
    class(globularClusterDissolutionRateDisksClass           ), intent(in   ), target :: globularClusterDissolutionRateDisks_
    class(globularClusterDissolutionRateSpheroidsClass       ), intent(in   ), target :: globularClusterDissolutionRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterDissolutionRateDisks_, *globularClusterDissolutionRateSpheroids_"/>
    !!]
    return
  end function globularClusterDissolutionRateConstructorInternal

  subroutine globularClusterDissolutionRateDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterDissolutionRate} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%globularClusterDissolutionRateDisks_"    />
    <objectDestructor name="self%globularClusterDissolutionRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterDissolutionRateDestructor

  double precision function globularClusterDissolutionRateExtract(self,node,time,instance)
    !!{
    Implement a star Dissolution rate output analysis property extractor.
    !!}
    implicit none
    double precision                                                     , dimension(:) , allocatable :: globularClusterDissolutionRateExtract
    class           (nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout), target      :: self
    type            (treeNode                                           ), intent(inout), target      :: node
    double precision                                                     , intent(in   )              :: time
    type            (multiCounter                                       ), intent(inout), optional    :: instance
    !$GLC attributes unused ::self, instance
    
    allocate(globularClusterDissolutionRateExtract(2))
    globularClusterDissolutionRateExtract(1)=+self%globularClusterDissolutionRateDisks_    %rate(node)
    globularClusterDissolutionRateExtract(2)=+self%globularClusterDissolutionRateSpheroids_%rate(node)
    return
  end function globularClusterDissolutionRateExtract

  subroutine globularClusterDissolutionRateNames(self,time,names)
    !!{
    Return the name of the globularClusterDissolutionRate property.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout)                             :: self
    double precision                                                     , intent(in   )                             :: time
    type            (varying_string                                     ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    allocate(names(2))
    names(1)=var_str('diskGlobularClusterDissolutionRate'    )
    names(2)=var_str('spheroidGlobularClusterDissolutionRate')
    return
  end subroutine globularClusterDissolutionRateNames

  subroutine globularClusterDissolutionRateDescriptions(self,time,descriptions)
    !!{
    Return a description of the globularClusterDissolutionRate property.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout)                             :: self
    double precision                                                     , intent(in   )                             :: time
    type            (varying_string                                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time
    allocate(descriptions(2))
    descriptions(1)=var_str('Disk globular cluster dissolution rate [M☉ Gyr⁻¹].')
    descriptions(2)=var_str('Spheroidal globular cluster dissolution rate [M☉ Gyr⁻¹].')
    return
  end subroutine globularClusterDissolutionRateDescriptions

  function globularClusterDissolutionRateUnitsInSI(self, time)
    !!{
    Return the units of the globularClusterDissolutionRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                                     , allocatable  , dimension(:) :: globularClusterDissolutionRateUnitsInSI
    class           (nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout)               :: self
    double precision                                                     , intent(in   )               :: time
    !$GLC attributes unused :: self, time
    allocate(globularClusterDissolutionRateUnitsInSI(2))
    globularClusterDissolutionRateUnitsInSI=[massSolar/gigaYear, massSolar/gigaYear]
    return
  end function globularClusterDissolutionRateUnitsInSI
