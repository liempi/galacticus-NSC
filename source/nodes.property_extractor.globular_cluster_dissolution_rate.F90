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
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorGlobularClusterDissolutionRate
     !!{
     A star Dissolution rate property extractor class.
     !!}
     private
     class(globularClusterDissolutionRateDisksClass    ), pointer :: globularClusterDissolutionRateDisks_    => null()
     class(globularClusterDissolutionRateSpheroidsClass), pointer :: globularClusterDissolutionRateSpheroids_=> null()
     type (varying_string                              )          :: name_                                            , description_, &
          &                                                        component
   contains
     final     ::                globularClusterDissolutionRateDestructor
     procedure :: extract     => globularClusterDissolutionRateExtract
     procedure :: name        => globularClusterDissolutionRateName
     procedure :: description => globularClusterDissolutionRateDescription
     procedure :: unitsInSI   => globularClusterDissolutionRateUnitsInSI
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
    type (varying_string                                     )                :: component
    type (enumerationGalacticComponentType                   )                :: component_

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star Dissolution rate.</description>
    </inputParameter>
    !!]
    component_=enumerationGalacticComponentEncode(char(component),includesPrefix=.false.)
    select case (component_%ID)
    case (galacticComponentDisk    %ID)
       !![
       <objectBuilder class="globularClusterDissolutionRateDisks"     name="globularClusterDissolutionRateDisks_"     source="parameters"/>
       !!]
       globularClusterDissolutionRateSpheroids_ => null()
    case (galacticComponentSpheroid%ID)
       !![
       <objectBuilder class="globularClusterDissolutionRateSpheroids" name="globularClusterDissolutionRateSpheroids_" source="parameters"/>
       !!]
       globularClusterDissolutionRateDisks_     => null()
    case (galacticComponentTotal   %ID)
       !![
       <objectBuilder class="globularClusterDissolutionRateDisks"     name="globularClusterDissolutionRateDisks_"     source="parameters"/>
       <objectBuilder class="globularClusterDissolutionRateSpheroids" name="globularClusterDissolutionRateSpheroids_" source="parameters"/>
       !!]
    end select
    self=nodePropertyExtractorGlobularClusterDissolutionRate(globularClusterDissolutionRateDisks_,globularClusterDissolutionRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterDissolutionRateDisks_"    />
    <objectDestructor name="globularClusterDissolutionRateSpheroids_"/>
    !!]
    return
  end function globularClusterDissolutionRateConstructorParameters

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

    if      (associated(self%globularClusterDissolutionRateDisks_).and.associated(self%globularClusterDissolutionRateSpheroids_)) then
       self%name_       ="totalglobularClusterDissolutionRate"
       self%description_="Total (disk + spheroid) globular cluster dissolution rate [M☉ Gyr⁻¹]."
       self%component   ="total"
    else if (associated(self%globularClusterDissolutionRateDisks_)                                                                                                            ) then
       self%name_       ="diskglobularClusterDissolutionRate"
       self%description_="Disk globular cluster dissolution rate [M☉ Gyr⁻¹]."
       self%component   ="disk"
    else if (                                                          associated(self%globularClusterDissolutionRateSpheroids_)                                                           ) then
       self%name_       ="spheroidglobularClusterDissolutionRate"
       self%description_="Spheroid globular cluster dissolution rate [M☉ Gyr⁻¹]."
       self%component   ="spheroid"
    else
       call Error_Report('No star Dissolution rate specified.'//{introspection:location})
    end if
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

  double precision function globularClusterDissolutionRateExtract(self,node,instance)
    !!{
    Implement a star Dissolution rate output analysis property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout), target   :: self
    type (treeNode                                           ), intent(inout), target   :: node
    type (multiCounter                                       ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    globularClusterDissolutionRateExtract=0.0d0
    if (associated(self%globularClusterDissolutionRateDisks_    )) globularClusterDissolutionRateExtract=globularClusterDissolutionRateExtract+self%globularClusterDissolutionRateDisks_    %rate(node)
    if (associated(self%globularClusterDissolutionRateSpheroids_)) globularClusterDissolutionRateExtract=globularClusterDissolutionRateExtract+self%globularClusterDissolutionRateSpheroids_%rate(node)
    return
  end function globularClusterDissolutionRateExtract

  function globularClusterDissolutionRateName(self)
    !!{
    Return the name of the globularClusterDissolutionRate property.
    !!}
    implicit none
    type (varying_string                                     )                :: globularClusterDissolutionRateName
    class(nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout) :: self

    globularClusterDissolutionRateName=self%name_
    return
  end function globularClusterDissolutionRateName

  function globularClusterDissolutionRateDescription(self)
    !!{
    Return a description of the globularClusterDissolutionRate property.
    !!}
    implicit none
    type (varying_string                                     )                :: globularClusterDissolutionRateDescription
    class(nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout) :: self

    globularClusterDissolutionRateDescription=self%description_
    return
  end function globularClusterDissolutionRateDescription

  double precision function globularClusterDissolutionRateUnitsInSI(self)
    !!{
    Return the units of the globularClusterDissolutionRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    class(nodePropertyExtractorGlobularClusterDissolutionRate), intent(inout) :: self
    !$GLC attributes unused :: self

    globularClusterDissolutionRateUnitsInSI=massSolar/gigaYear
    return
  end function globularClusterDissolutionRateUnitsInSI
