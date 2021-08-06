!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a property extractor class for the orbital adiabatic ratio of disks.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorAdiabaticRatioOrbitalDisk">
    <description>
      A property extractor class for the orbital adiabatic ratio of disks. The orbital adiabatic ratio, $\mathcal{R}$, is defined
      as:
      \begin{equation}
      \mathcal{R} = \frac{r_\mathrm{p} / v_\mathrm{p}}{2 \pi R_\mathrm{d} / v_\mathrm{d}},
      \end{equation}      
      where $r_\mathrm{p}$ and $v_\mathrm{p}$ are the orbital radius and velocity at pericenter respectively, and $R_\mathrm{d}$
      and $v_\mathrm{d}$ are the characteristic radius of the disk and the rotation curve at that radius respectively.
    </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorAdiabaticRatioOrbitalDisk
     !!{
     A property extractor class for the orbital adiabatic ratio of disks.
     !!}
     private
   contains
     procedure :: extract     => adiabaticRatioOrbitalDiskExtract
     procedure :: name        => adiabaticRatioOrbitalDiskName
     procedure :: description => adiabaticRatioOrbitalDiskDescription
     procedure :: unitsInSI   => adiabaticRatioOrbitalDiskUnitsInSI
     procedure :: type        => adiabaticRatioOrbitalDiskType
  end type nodePropertyExtractorAdiabaticRatioOrbitalDisk

  interface nodePropertyExtractorAdiabaticRatioOrbitalDisk
     !!{
     Constructors for the ``adiabaticRatioOrbitalDisk'' output analysis class.
     !!}
     module procedure adiabaticRatioOrbitalDiskConstructorParameters
  end interface nodePropertyExtractorAdiabaticRatioOrbitalDisk

contains

  function adiabaticRatioOrbitalDiskConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily adiabaticRatioOrbitalDisk} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorAdiabaticRatioOrbitalDisk)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    self=nodePropertyExtractorAdiabaticRatioOrbitalDisk()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function adiabaticRatioOrbitalDiskConstructorParameters

  double precision function adiabaticRatioOrbitalDiskExtract(self,node,instance)
    !!{
    Implement extraction of the orbital adiabatic ratio for disks.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentDisk                               , nodeComponentSatellite
    use :: Kepler_Orbits           , only : keplerOrbit
    use :: Satellite_Orbits        , only : Satellite_Orbit_Extremum_Phase_Space_Coordinates, extremumPericenter
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (nodePropertyExtractorAdiabaticRatioOrbitalDisk), intent(inout)           :: self
    type            (treeNode                                      ), intent(inout), target   :: node
    type            (multiCounter                                  ), intent(inout), optional :: instance
    class           (nodeComponentDisk                             )               , pointer  :: disk
    class           (nodeComponentSatellite                        )               , pointer  :: satellite
    type            (treeNode                                      )               , pointer  :: nodeHost
    type            (keplerOrbit                                   )                          :: orbit
    double precision                                                                          :: radiusPericenter, velocityPericenter
    !$GLC attributes unused :: instance

    adiabaticRatioOrbitalDiskExtract=-1.0d0
    if (.not.node%isSatellite()) return
    nodeHost  => node     %parent
    disk      => node     %disk       ()
    satellite => node     %satellite  ()
    orbit     =  satellite%virialOrbit()
    call Satellite_Orbit_Extremum_Phase_Space_Coordinates(nodeHost,orbit,extremumPericenter,radiusPericenter,velocityPericenter)
    if (disk%radius() > 0.0d0)                                    &
         & adiabaticRatioOrbitalDiskExtract=+(                    &
         &                                    +  radiusPericenter &
         &                                    /velocityPericenter &
         &                                   )                    &
         &                                  /(                    &
         &                                    +2.0d0              &
         &                                    *Pi                 &
         &                                    *disk%radius  ()    &
         &                                    /disk%velocity()    &
         &                                   )
    return
  end function adiabaticRatioOrbitalDiskExtract

  function adiabaticRatioOrbitalDiskName(self)
    !!{
    Return the name of the orbital adiabatic ratio for disks property.
    !!}
    implicit none
    type (varying_string                                )                :: adiabaticRatioOrbitalDiskName
    class(nodePropertyExtractorAdiabaticRatioOrbitalDisk), intent(inout) :: self
    !$GLC attributes unused :: self

    adiabaticRatioOrbitalDiskName=var_str('ratioAdiabaticOrbitalDisk')
    return
  end function adiabaticRatioOrbitalDiskName

  function adiabaticRatioOrbitalDiskDescription(self)
    !!{
    Return a description of the orbital adiabatic ratio for disks property.
    !!}
    implicit none
    type (varying_string                                )                :: adiabaticRatioOrbitalDiskDescription
    class(nodePropertyExtractorAdiabaticRatioOrbitalDisk), intent(inout) :: self
    !$GLC attributes unused :: self

    adiabaticRatioOrbitalDiskDescription=var_str('The orbital adiabatic ratio of the disk.')
    return
  end function adiabaticRatioOrbitalDiskDescription

  double precision function adiabaticRatioOrbitalDiskUnitsInSI(self)
    !!{
    Return the units of the orbital adiabatic ratio for disks property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorAdiabaticRatioOrbitalDisk), intent(inout) :: self
    !$GLC attributes unused :: self

    adiabaticRatioOrbitalDiskUnitsInSI=1.0d0
    return
  end function adiabaticRatioOrbitalDiskUnitsInSI

  integer function adiabaticRatioOrbitalDiskType(self)
    !!{
    Return the type of the orbital adiabatic ratio for disks property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorAdiabaticRatioOrbitalDisk), intent(inout) :: self
    !$GLC attributes unused :: self

    adiabaticRatioOrbitalDiskType=outputAnalysisPropertyTypeLinear
    return
  end function adiabaticRatioOrbitalDiskType
