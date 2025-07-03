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
Implements a globular cluster dissolution rate property extractor class.
!!}

  use :: Globular_Cluster_Infall_Rates_Disks    , only : globularClusterInfallRateDisksClass
  use :: Globular_Cluster_Infall_Rates_Spheroids, only : globularClusterInfallRateSpheroidsClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorGlobularClusterInfallRate">
   <description>
    A node property extractor which extracts the globular cluster infall rate in a galaxy. The type of globular dissolution rate is controlled by
    the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'' or ``{\normalfont \ttfamily total}''. The corresponding globular cluster infall
    rate is extracted as {\normalfont \ttfamily \textless\ component\textgreater\ globularClusterInfallRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorGlobularClusterInfallRate
     !!{
     A star Dissolution rate property extractor class.
     !!}
     private
     class(globularClusterInfallRateDisksClass    ), pointer :: globularClusterInfallRateDisks_    => null()
     class(globularClusterInfallRateSpheroidsClass), pointer :: globularClusterInfallRateSpheroids_=> null()
   contains
     final     ::                 globularClusterInfallRateDestructor
     procedure :: elementCount => globularClusterInfallRateElementCount
     procedure :: extract      => globularClusterInfallRateExtract
     procedure :: names        => globularClusterInfallRateNames
     procedure :: descriptions => globularClusterInfallRateDescriptions
     procedure :: unitsInSI    => globularClusterInfallRateUnitsInSI
  end type nodePropertyExtractorGlobularClusterInfallRate

  interface nodePropertyExtractorGlobularClusterInfallRate
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterInfallRate} output analysis class.
     !!}
     module procedure globularClusterInfallRateConstructorParameters
     module procedure globularClusterInfallRateConstructorInternal
  end interface nodePropertyExtractorGlobularClusterInfallRate

contains

  function globularClusterInfallRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterInfallRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorGlobularClusterInfallRate)                :: self
    type (inputParameters                               ), intent(inout) :: parameters
    class(globularClusterInfallRateDisksClass           ), pointer       :: globularClusterInfallRateDisks_
    class(globularClusterInfallRateSpheroidsClass       ), pointer       :: globularClusterInfallRateSpheroids_

    !![
    <objectBuilder class="globularClusterInfallRateDisks"     name="globularClusterInfallRateDisks_"     source="parameters"/>
    <objectBuilder class="globularClusterInfallRateSpheroids" name="globularClusterInfallRateSpheroids_" source="parameters"/>
    !!]
    self=nodePropertyExtractorGlobularClusterInfallRate(globularClusterInfallRateDisks_,globularClusterInfallRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterInfallRateDisks_"    />
    <objectDestructor name="globularClusterInfallRateSpheroids_"/>
    !!]
    return
  end function globularClusterInfallRateConstructorParameters

  integer function globularClusterInfallRateElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily radiiHalfLightProperties} property extractor class.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterInfallRate), intent(inout) :: self
    double precision                                                , intent(in   ) :: time
    !$GLC attributes unused :: self

    globularClusterInfallRateElementCount=2
    return
  end function globularClusterInfallRateElementCount

  function globularClusterInfallRateConstructorInternal(globularClusterInfallRateDisks_,globularClusterInfallRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterInfallRate} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorGlobularClusterInfallRate)                        :: self
    class(globularClusterInfallRateDisksClass           ), intent(in   ), target :: globularClusterInfallRateDisks_
    class(globularClusterInfallRateSpheroidsClass       ), intent(in   ), target :: globularClusterInfallRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterInfallRateDisks_, *globularClusterInfallRateSpheroids_"/>
    !!]
    return
  end function globularClusterInfallRateConstructorInternal

  subroutine globularClusterInfallRateDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterInfallRate} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGlobularClusterInfallRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%globularClusterInfallRateDisks_"    />
    <objectDestructor name="self%globularClusterInfallRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterInfallRateDestructor

  double precision function globularClusterInfallRateExtract(self,node,time,instance)
    !!{
    Implement a star Dissolution rate output analysis property extractor.
    !!}
    implicit none
    double precision                                                , dimension(:) , allocatable :: globularClusterInfallRateExtract
    class           (nodePropertyExtractorGlobularClusterInfallRate), intent(inout), target      :: self
    type            (treeNode                                      ), intent(inout), target      :: node
    double precision                                                , intent(in   )              :: time
    type            (multiCounter                                  ), intent(inout), optional    :: instance
    !$GLC attributes unused ::self, instance
    allocate(globularClusterInfallRateExtract(2))
    globularClusterInfallRateExtract(1)=+self%globularClusterInfallRateDisks_    %rate(node)
    globularClusterInfallRateExtract(2)=+self%globularClusterInfallRateSpheroids_%rate(node)
    return
  end function globularClusterInfallRateExtract

  subroutine globularClusterInfallRateNames(self,time,names)
    !!{
    Return the name of the globularClusterInfallRate property.
    !!}
    implicit none
    class           (nodePropertyExtractorGlobularClusterInfallRate), intent(inout)                             :: self
    double precision                                                , intent(in   )                             :: time
    type            (varying_string                                ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    allocate(names(2))
    names(1)=var_str('diskGlobularClusterInfallRate'    )
    names(2)=var_str('spheroidGlobularClusterInfallRate')
    return
  end subroutine globularClusterInfallRateNames

  subroutine globularClusterInfallRateDescriptions(self,time,descriptions)
    !!{
    Return a description of the globularClusterInfallRate property.
    !!}
    implicit none
    class(nodePropertyExtractorGlobularClusterInfallRate), intent(inout) :: self
    double precision                                     , intent(in   )                             :: time
    type            (varying_string                     ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time
    allocate(descriptions(2))
    descriptions(1)=var_str('Disk globular clusters infall rate [M☉ Gyr⁻¹].')
    descriptions(2)=var_str('Spheroidal globular clusters infall rate [M☉ Gyr⁻¹].')
    return
  end subroutine globularClusterInfallRateDescriptions

  double precision function globularClusterInfallRateUnitsInSI(self, time)
    !!{
    Return the units of the globularClusterInfallRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                                , allocatable  , dimension(:) :: globularClusterInfallRateUnitsInSI
    class           (nodePropertyExtractorGlobularClusterInfallRate), intent(inout)               :: self
    double precision                                                , intent(in   )               :: time
    !$GLC attributes unused :: self, time
    allocate(globularClusterInfallRateUnitsInSI(2))
    globularClusterInfallRateUnitsInSI=[massSolar/gigaYear, massSolar/gigaYear]
    return
  end function globularClusterInfallRateUnitsInSI
