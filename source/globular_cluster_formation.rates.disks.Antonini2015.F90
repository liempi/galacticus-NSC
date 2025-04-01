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
  Implementation of a star formation rate in galactic disks which computes the rate by integrating a star formation rate over
  the disk.
  !!}
  use :: Star_Formation_Rates_Disks, only : starFormationRateDisksClass

  !![
  <globularClusterFormationRateDisks name="globularClusterFormationRateDisksAntonini2015">
   <description>
    A globular cluster formation rate in galactic disks which computes the rate by multiply the star formation rate of the disk by a factor. Specifically, the globular cluster formation rate is given by
    \begin{equation}
     \dot{M}_\mathrm{disk}^\mathrm{gc} = f_\mathrm{gc} \dot{\Sigma}_\star(r) \mathrm{d}r,
    \end{equation}
    where $\dot{\Sigma}_\star(r)$ is the surface density of star formation rate.
   </description>
  </globularClusterFormationRateDisks>
  !!]
  type, extends(globularClusterFormationRateDisksClass) :: globularClusterFormationRateDisksAntonini2015
     !!{
     Implementation of a rate for globular cluster formation in galactic disks which computes the rate by multiplyiing the star formation rate
     of the disk by an efficiency parameter.
     !!}
     private
     class           (starFormationRateDisksClass), pointer :: starFormationRateDisks_ => null()
     double precision                                       :: efficiency                              
   contains
     final     ::         globularClusterFormationDisksAntonini2015Destructor
     procedure :: rate => globularClusterFormationDisksAntonini2015Rate
  end type globularClusterFormationRateDisksAntonini2015

  interface globularClusterFormationRateDisksAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationDisksAntonini2015} formation rate in disks class.
     !!}
     module procedure globularClusterFormationDisksAntonini2015ConstructorParameters
     module procedure globularClusterFormationDisksAntonini2015ConstructorInternal
  end interface globularClusterFormationRateDisksAntonini2015
    
contains

  function globularClusterFormationDisksAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormationDisks} formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterFormationRateDisksAntonini2015)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (starFormationRateDisksClass                  ), pointer       :: starFormationRateDisks_
    double precision                                                               :: efficiency

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultValue>7.0d-2</defaultValue>
      <description>Relative tolerance to use when integrating star formation rate surface densities over the disk.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateDisks" name="starFormationRateDisks_" source="parameters"/>
    !!]
    self=globularClusterFormationRateDisksAntonini2015(efficiency,starFormationRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"/>
    !!]
    return
  end function globularClusterFormationDisksAntonini2015ConstructorParameters

  function globularClusterFormationDisksAntonini2015ConstructorInternal(efficiency,starFormationRateDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationDisks} globular cluster formation rate in disks class.
    !!}
    implicit none
    type            (globularClusterFormationRateDisksAntonini2015)                        :: self
    class           (starFormationRateDisksClass                  ), intent(in   ), target :: starFormationRateDisks_
    double precision                                               , intent(in   )         :: efficiency
    !![
    <constructorAssign variables="efficiency,*starFormationRateDisks_"/>
    !!]
    return
  end function globularClusterFormationDisksAntonini2015ConstructorInternal

  subroutine globularClusterFormationDisksAntonini2015Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationDisks} class
    !!}
    implicit none
    type(globularClusterFormationRateDisksAntonini2015), intent(inout) :: self
    !![
    <objectDestructor name="self%starFormationRateDisks_"/>
    !!]
    return
  end subroutine globularClusterFormationDisksAntonini2015Destructor

  double precision function globularClusterFormationDisksAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily node}
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentDisk, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (globularClusterFormationRateDisksAntonini2015), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    class           (nodeComponentDisk                            ), pointer                :: disk
    double precision                                                                        :: radiusDisk           , massStellar, &
         &                                                                                     rateStarFormationDisk

    ! Get the disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    radiusDisk  =  disk%radius     ()
    ! Check if the disk is physical.
    if (massStellar <= 0.0d0 .or. radiusDisk <= 0.0d0) then
        ! It is not, so return zero rate.
        rate=+0.0d0
    else
        rateStarFormationDisk =  self%starFormationRateDisks_%rate    (node)  
        ! Find the rate of globular cluster formation in the disk component.
        if (rateStarFormationDisk <= 0.0d0) then
            rate   =+0.0d0
        else
            ! Gas accretion rate model from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
            rate   =+self%efficiency                &
                 &  *     rateStarFormationDisk
        end if 
    end if
    return
  end function globularClusterFormationDisksAntonini2015Rate

  
