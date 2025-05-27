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
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorGlobularClusterInfallRate
     !!{
     A star Dissolution rate property extractor class.
     !!}
     private
     class(globularClusterInfallRateDisksClass    ), pointer :: globularClusterInfallRateDisks_    => null()
     class(globularClusterInfallRateSpheroidsClass), pointer :: globularClusterInfallRateSpheroids_=> null()
     type (varying_string                         )          :: name_                                            , description_, &
          &                                                    component
   contains
     final     ::                globularClusterInfallRateDestructor
     procedure :: extract     => globularClusterInfallRateExtract
     procedure :: name        => globularClusterInfallRateName
     procedure :: description => globularClusterInfallRateDescription
     procedure :: unitsInSI   => globularClusterInfallRateUnitsInSI
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
    type (varying_string                                )                :: component
    type (enumerationGalacticComponentType              )                :: component_

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
       <objectBuilder class="globularClusterInfallRateDisks"     name="globularClusterInfallRateDisks_"     source="parameters"/>
       !!]
       globularClusterInfallRateSpheroids_ => null()
    case (galacticComponentSpheroid%ID)
       !![
       <objectBuilder class="globularClusterInfallRateSpheroids" name="globularClusterInfallRateSpheroids_" source="parameters"/>
       !!]
       globularClusterInfallRateDisks_     => null()
    case (galacticComponentTotal   %ID)
       !![
       <objectBuilder class="globularClusterInfallRateDisks"     name="globularClusterInfallRateDisks_"     source="parameters"/>
       <objectBuilder class="globularClusterInfallRateSpheroids" name="globularClusterInfallRateSpheroids_" source="parameters"/>
       !!]
    end select
    self=nodePropertyExtractorGlobularClusterInfallRate(globularClusterInfallRateDisks_,globularClusterInfallRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterInfallRateDisks_"    />
    <objectDestructor name="globularClusterInfallRateSpheroids_"/>
    !!]
    return
  end function globularClusterInfallRateConstructorParameters

  function globularClusterInfallRateConstructorInternal(globularClusterInfallRateDisks_,globularClusterInfallRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterInfallRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorGlobularClusterInfallRate)                        :: self
    class(globularClusterInfallRateDisksClass           ), intent(in   ), target :: globularClusterInfallRateDisks_
    class(globularClusterInfallRateSpheroidsClass       ), intent(in   ), target :: globularClusterInfallRateSpheroids_
    !![
    <constructorAssign variables="*globularClusterInfallRateDisks_, *globularClusterInfallRateSpheroids_"/>
    !!]

    if      (associated(self%globularClusterInfallRateDisks_).and.associated(self%globularClusterInfallRateSpheroids_)) then
       self%name_       ="totalglobularClusterInfallRate"
       self%description_="Total (disk + spheroid) globular cluster infall rate [M☉ Gyr⁻¹]."
       self%component   ="total"
    else if (associated(self%globularClusterInfallRateDisks_)                                                                                                            ) then
       self%name_       ="diskglobularClusterInfallRate"
       self%description_="Disk globular cluster infall rate [M☉ Gyr⁻¹]."
       self%component   ="disk"
    else if (                                                          associated(self%globularClusterInfallRateSpheroids_)                                                           ) then
       self%name_       ="spheroidglobularClusterInfallRate"
       self%description_="Spheroid globular cluster infall rate [M☉ Gyr⁻¹]."
       self%component   ="spheroid"
    else
       call Error_Report('No star infall rate specified.'//{introspection:location})
    end if
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

  double precision function globularClusterInfallRateExtract(self,node,instance)
    !!{
    Implement a star Dissolution rate output analysis property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorGlobularClusterInfallRate), intent(inout), target   :: self
    type (treeNode                                           ), intent(inout), target   :: node
    type (multiCounter                                       ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    globularClusterInfallRateExtract=0.0d0
    if (associated(self%globularClusterInfallRateDisks_    )) globularClusterInfallRateExtract=globularClusterInfallRateExtract+self%globularClusterInfallRateDisks_    %rate(node)
    if (associated(self%globularClusterInfallRateSpheroids_)) globularClusterInfallRateExtract=globularClusterInfallRateExtract+self%globularClusterInfallRateSpheroids_%rate(node)
    return
  end function globularClusterInfallRateExtract

  function globularClusterInfallRateName(self)
    !!{
    Return the name of the globularClusterInfallRate property.
    !!}
    implicit none
    type (varying_string                                )                :: globularClusterInfallRateName
    class(nodePropertyExtractorGlobularClusterInfallRate), intent(inout) :: self

    globularClusterInfallRateName=self%name_
    return
  end function globularClusterInfallRateName

  function globularClusterInfallRateDescription(self)
    !!{
    Return a description of the globularClusterInfallRate property.
    !!}
    implicit none
    type (varying_string                                )                :: globularClusterInfallRateDescription
    class(nodePropertyExtractorGlobularClusterInfallRate), intent(inout) :: self

    globularClusterInfallRateDescription=self%description_
    return
  end function globularClusterInfallRateDescription

  double precision function globularClusterInfallRateUnitsInSI(self)
    !!{
    Return the units of the globularClusterInfallRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    class(nodePropertyExtractorGlobularClusterInfallRate), intent(inout) :: self
    !$GLC attributes unused :: self

    globularClusterInfallRateUnitsInSI=massSolar/gigaYear
    return
  end function globularClusterInfallRateUnitsInSI
