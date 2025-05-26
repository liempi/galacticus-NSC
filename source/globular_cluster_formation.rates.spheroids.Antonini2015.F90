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
  Implementation of a star formation rate in galactic spheroids which computes the rate by integrating a star formation rate over
  the spheroid.
  !!}
  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass

  !![
  <globularClusterFormationRateSpheroids name="globularClusterFormationRateSpheroidsAntonini2015">
   <description>
    A globular cluster formation rate in galactic spheroids which computes the rate by multiply the star formation rate of the spheroid by a factor. Specifically, the globular cluster formation rate in
    spheroids is given by
    \begin{equation}
     \dot{M}_\mathrm{spheroid}^\mathrm{gc} = f_\mathrm{gc} \dot{M}_\star(r),
    \end{equation}
    where $\dot{M}_\star$ is the star formation rate of the spheroid.
   </description>
  </globularClusterFormationRateSpheroids>
  !!]
  type, extends(globularClusterFormationRateSpheroidsClass) :: globularClusterFormationRateSpheroidsAntonini2015
     !!{
     Implementation of a rate for globular cluster formation in galactic spheroids which computes the rate by multiplyiing the star formation rate
     of the spheroid by an efficiency parameter.
     !!}
     private
     class           (starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
     double precision                                           :: efficiency                              
   contains
     final     ::         globularClustersSpheroidsAntonini2015Destructor
     procedure :: rate => globularClustersSpheroidsAntonini2015Rate
  end type globularClusterFormationRateSpheroidsAntonini2015

  interface globularClusterFormationRateSpheroidsAntonini2015
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationSpheroids} formation rate in spheroids class.
     !!}
     module procedure globularClustersSpheroidsAntonini2015ConstructorParameters
     module procedure globularClustersSpheroidsAntonini2015ConstructorInternal
  end interface globularClusterFormationRateSpheroidsAntonini2015
    
contains

  function globularClustersSpheroidsAntonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} formation rate in spheroids class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (globularClusterFormationRateSpheroidsAntonini2015)                :: self
    type            (inputParameters                                  ), intent(inout) :: parameters
    class           (starFormationRateSpheroidsClass                  ), pointer       :: starFormationRateSpheroids_
    double precision                                                                   :: efficiency

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultValue>7.0d-2</defaultValue>
      <description>Efficiency of the globular cluster formation rate in the spheroid component.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=globularClusterFormationRateSpheroidsAntonini2015(efficiency,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function globularClustersSpheroidsAntonini2015ConstructorParameters

  function globularClustersSpheroidsAntonini2015ConstructorInternal(efficiency,starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} globular cluster formation rate in spheroids class.
    !!}
    implicit none
    type            (globularClusterFormationRateSpheroidsAntonini2015)                        :: self
    class           (starFormationRateSpheroidsClass                  ), intent(in   ), target :: starFormationRateSpheroids_
    double precision                                                   , intent(in   )         :: efficiency
    !![
    <constructorAssign variables="efficiency,*starFormationRateSpheroids_"/>
    !!]
    return
  end function globularClustersSpheroidsAntonini2015ConstructorInternal

  subroutine globularClustersSpheroidsAntonini2015Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} class
    !!}
    implicit none
    type(globularClusterFormationRateSpheroidsAntonini2015), intent(inout) :: self
    !![
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    !!]
    return
  end subroutine globularClustersSpheroidsAntonini2015Destructor

  double precision function globularClustersSpheroidsAntonini2015Rate(self,node) result(rate)
    !!{
    Returns the globular formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic spheroid of {\normalfont \ttfamily node}
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid, treeNode
    implicit none
    class           (globularClusterFormationRateSpheroidsAntonini2015), intent(inout), target  :: self
    type            (treeNode                                         ), intent(inout), target  :: node
    class           (nodeComponentSpheroid                            ), pointer                :: spheroid
    double precision                                                                            :: radiusSpheroid           , massStellar, &
         &                                                                                         rateStarFormationSpheroid

    ! Get the spheroid properties.
    spheroid       => node    %spheroid   ()
    massStellar    =  spheroid%massStellar()
    radiusSpheroid =  spheroid%radius     ()

    ! Check if the spheroid is physical.
    if (massStellar <= 0.0d0 .or. radiusSpheroid <= 0.0d0) then
        ! It is not, so return zero rate.
        rate=+0.0d0
    else
        rateStarFormationSpheroid =  self%starFormationRateSpheroids_%rate    (node)  
        ! Find the rate of globular cluster formation in the spheroid component.
        if (rateStarFormationSpheroid <= 0.0d0) then
            rate   =+0.0d0
        else
            ! Globular cluster formation rate from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
            rate   =+self%efficiency                &
                 &  *     rateStarFormationSpheroid
        end if 
    end if
    return
  end function globularClustersSpheroidsAntonini2015Rate
