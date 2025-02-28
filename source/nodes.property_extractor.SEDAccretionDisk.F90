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
  Contains a module which implements a property extractor class for the SED of the accreting disk.
  !!}
  use :: Cosmology_Functions   , only : cosmologyFunctionsClass
  use :: Accretion_Disk_Spectra, only : accretionDiskSpectraClass
     
  !![
  <nodePropertyExtractor name="nodePropertyExtractorSEDaccretionDisk">
    <description>A property extractor class for the SED of the accreting disk.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorSEDaccretionDisk
     !!{
     A property extractor class for the SED of a component.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer                   :: cosmologyFunctions_   => null()
     class           (accretionDiskSpectraClass), pointer                   :: accretionDiskSpectra_ => null()
     integer         (c_size_t                 )                            :: countWavelengths=-1_c_size_t
     double precision                           , allocatable, dimension(:) :: wavelengths_                                    
     double precision                                                       :: wavelengthMinimum               , wavelengthMaximum, &
          &                                                                    resolution                      , factorWavelength 
     type            (enumerationFrameType     )                            :: frame
   contains
     !![
     <methods>
       <method description="Return an array of the wavelengths at which the SED is computed."  method="wavelengths"/>
     </methods>
     !!]
     final     ::                            sedAccretionDiskDestructor
     procedure :: wavelengths             => sedAccretionDiskWavelengths
     procedure :: columnDescriptions      => sedAccretionDiskColumnDescriptions
     procedure :: size                    => sedAccretionDiskSize
     procedure :: elementCount            => sedAccretionDiskElementCount
     procedure :: extract                 => sedAccretionDiskExtract
     procedure :: names                   => sedAccretionDiskNames
     procedure :: descriptions            => sedAccretionDiskDescriptions
     procedure :: unitsInSI               => sedAccretionDiskUnitsInSI
  end type nodePropertyExtractorSEDaccretionDisk
  
  interface nodePropertyExtractorSEDaccretionDisk
     !!{
     Constructors for the ``sed'' output analysis class.
     !!}
     module procedure sedAccretionDiskConstructorParameters
     module procedure sedAccretionDiskConstructorInternal
  end interface nodePropertyExtractorSEDaccretionDisk
      
contains

  function sedAccretionDiskConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sed} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters              , only : inputParameter, inputParameters
    use :: Stellar_Luminosities_Structure, only : enumerationFrameEncode

    implicit none
    type            (nodePropertyExtractorSEDaccretionDisk)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (accretionDiskSpectraClass            ), pointer       :: accretionDiskSpectra_
    type            (varying_string                       )                :: frame
    double precision                                                       :: wavelengthMinimum    , wavelengthMaximum, &
         &                                                                    resolution                            
    
    !![
    <inputParameter>
      <name>frame</name>
      <source>parameters</source>
      <defaultValue>var_str('rest')</defaultValue>
      <description>The frame ({\normalfont \ttfamily rest} or {\normalfont \ttfamily observed}) for which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum wavelength at which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum wavelength at which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>resolution</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>The resolution, $\lambda/\Delta\lambda$, at which to compute the SED. If a negative value is given the SED will be computed at the full resolution provided by the stellar population spectra class.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="accretionDiskSpectra"  name="accretionDiskSpectra_" source="parameters"/>
    !!]
    self=nodePropertyExtractorSEDaccretionDisk(enumerationFrameEncode(char(frame),includesPrefix=.false.),wavelengthMinimum,wavelengthMaximum,resolution,cosmologyFunctions_,accretionDiskSpectra_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="accretionDiskSpectra_"/>
    !!]
    return
  end function sedAccretionDiskConstructorParameters

  function sedAccretionDiskConstructorInternal(frame,wavelengthMinimum,wavelengthMaximum,resolution,cosmologyFunctions_,accretionDiskSpectra_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorSEDaccretionDisk)                              :: self
    type            (enumerationFrameType                 ), intent(in   )               :: frame
    class           (cosmologyFunctionsClass              ), intent(in   ), target       :: cosmologyFunctions_
    class           (accretionDiskSpectraClass            ), intent(in   ), target       :: accretionDiskSpectra_
    double precision                                       , intent(in   )               :: wavelengthMinimum    , wavelengthMaximum, &
         &                                                                                  resolution                            
    !![
    <constructorAssign variables="frame,wavelengthMinimum,wavelengthMaximum,resolution,*cosmologyFunctions_,*accretionDiskSpectra_"/>
    !!]
    
    ! Compute the factor by which the minimum/maximum wavelength in a resolution element differ from the central wavelength.
    if (resolution > 0.0d0) self%factorWavelength=(1.0d0+sqrt(1.0d0+4.0d0*resolution**2))/2.0d0/resolution
    return
  end function sedAccretionDiskConstructorInternal

  subroutine sedAccretionDiskDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSEDaccretionDisk), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%accretionDiskSpectra_"/>
    !!]
    return
  end subroutine sedAccretionDiskDestructor

  integer function sedAccretionDiskElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily sed} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    sedAccretionDiskElementCount=1
    return
  end function sedAccretionDiskElementCount

  function sedAccretionDiskSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily sed} property extractors.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved 
    implicit none
    integer         (c_size_t                             )                              :: sedAccretionDiskSize
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)               :: self
    double precision                                       , intent(in   )               :: time
    logical                                                , allocatable  , dimension(:) :: selection
    double precision                                                                     :: expansionFactor

    if (self%resolution < 0.0d0) then
       ! Full resolution, so the number of wavelengths is simply the total number available within the wavelength range.
       allocate(selection(self%countWavelengths))
       select case (self%frame%ID)
       case (frameRest    %ID)
          expansionFactor=1.0d0
          selection      =                                                             &
               &           self%wavelengths_                 >= self%wavelengthMinimum &
               &          .and.                                                        &
               &           self%wavelengths_                 <= self%wavelengthMaximum
       case (frameObserved%ID)
          expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
          selection      =                                                             &
               &           self%wavelengths_/expansionFactor >= self%wavelengthMinimum &
               &          .and.                                                        &
               &           self%wavelengths_/expansionFactor <= self%wavelengthMaximum
       case default
          expansionFactor      =1.0d0
          sedAccretionDiskSize =0_c_size_t
          call Error_Report('unknown frame'//{introspection:location})
       end select
       sedAccretionDiskSize      =count(selection)
    else
       ! Compute the number of wavelengths.
       sedAccretionDiskSize=int(log(self%wavelengthMaximum/self%wavelengthMinimum)/log(self%factorWavelength)/2.0d0)+1
    end if
    return
  end function sedAccretionDiskSize

  function sedAccretionDiskExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily sed} property extractor.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    double precision                                       , dimension(:,:  )          , allocatable :: sedAccretionDiskExtract
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)   , target                :: self
    type            (treeNode                             ), intent(inout)   , target                :: node
    double precision                                       , intent(in   )                           :: time
    type            (multiCounter                         ), intent(inout)   , optional              :: instance
    integer                                                                                          :: iWavelength
    !$GLC attributes unused :: instance

    allocate(sedAccretionDiskExtract(self%size(time),1))

    sedAccretionDiskExtract=0.0d0
    do iWavelength=1,size(sedAccretionDiskExtract,dim=1)
       sedAccretionDiskExtract(iWavelength,1)=self%accretionDiskSpectra_%spectrumNode(node,self%wavelengths_(iWavelength))
    end do
    return
  end function sedAccretionDiskExtract

  subroutine sedAccretionDiskNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily sed} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(1))
    names(1)="AccretionDiskSED:"
    return
  end subroutine sedAccretionDiskNames

  subroutine sedAccretionDiskDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily sed} property.
    !!}
    implicit none
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(1))
    descriptions(1)="Spectral energy density (SED), dL/dν for the accreting disk [L☉ Hz⁻¹]."
    return
  end subroutine sedAccretionDiskDescriptions

  function sedAccretionDiskWavelengths(self,time)
    !!{
    Return wavelengths at which the SED is tabulated.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved
    implicit none
    double precision                                       , dimension(:) , allocatable :: sedAccretionDiskWavelengths
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    integer         (c_size_t                             )                             :: i             , j              
    double precision                                                                    :: wavelength    , expansionFactor

    allocate(sedAccretionDiskWavelengths(self%size(time)))
    j =0
    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=0.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select

    do i=1,size(sedAccretionDiskWavelengths)
       if (self%resolution < 0.0d0) then
          ! Full resolution SED.
          j=j+1
          do while (self%wavelengths_(j)/expansionFactor < self%wavelengthMinimum)
             j=j+1
          end do
          wavelength=self%wavelengths_(j)/expansionFactor
       else
          ! Finite resolution SED.
          if (i == 1) then
             wavelength=+self%wavelengthMinimum    &
                  &     *self%factorWavelength
          else
             wavelength=+     wavelength           &
                  &     *self%factorWavelength **2
          end if
       end if
       sedAccretionDiskWavelengths(i)=wavelength
    end do
    return
  end function sedAccretionDiskWavelengths

  subroutine sedAccretionDiskColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily sed} property.
    !!}
    use :: Numerical_Constants_Units, only : metersToAngstroms
    implicit none
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)                            :: self
    double precision                                       , intent(in   ), optional                  :: time
    type            (varying_string                       ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                       , intent(inout), dimension(:), allocatable :: values 
    type            (varying_string                       ), intent(  out)                            :: valuesDescription
    double precision                                       , intent(  out)                            :: valuesUnitsInSI
    integer         (c_size_t                             )                                           :: i
    character       (len=18                               )                                           :: label
    
    allocate(descriptions(self%size(time)))
    allocate(values      (self%size(time)))
    values=self%wavelengths(time)
    do i=1,size(descriptions)      
       write (label,'(a2,1x,e12.6,1x,a1)') "λ=",values(i),"Å"
       descriptions(i)=trim(label)
    end do
    valuesDescription=var_str('Wavelengths at which the SED is tabulated [in units of Å].')
    valuesUnitsInSI  =1.0d0/metersToAngstroms
    return
  end subroutine sedAccretionDiskColumnDescriptions

  function sedAccretionDiskUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily sed} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    double precision                                       , allocatable  , dimension(:) :: sedAccretionDiskUnitsInSI
    class           (nodePropertyExtractorSEDaccretionDisk), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(sedAccretionDiskUnitsInSI(1))
    sedAccretionDiskUnitsInSI(1)=luminositySolar
    return
  end function sedAccretionDiskUnitsInSI 