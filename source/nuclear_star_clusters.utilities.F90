!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Contains a module that implements various useful utility functions for calculations for the nuclear star cluster component.
  !!}

module Nuclear_Star_Clusters_Utilities
  !!{
  Implements various useful utility functions for calculations for the nuclear star cluster component.
  !!}
  implicit none
  private
  public :: nuclearStarClusterCrossingTimescale    , nuclearStarClusterKroupaNumberOfStars, &
    &       nuclearStarClusterCoreCollapseTimescale, nuclearStarClusterRelaxationTimescale

  contains

  double precision function nuclearStarClusterKroupaNumberOfStars(massStellar) result(numberOfStars) 
    !!{
       Returns the number of stars assuming a Kroupa IMF, assuming a minimum star mass M0.001, and maximum M125.
    !!}
    implicit none
    double precision, intent(in   ) :: massStellar

    numberOfStars = +0.380d0     &
       &            *massStellar &
       &            /0.079d0
    return
  end function nuclearStarClusterKroupaNumberOfStars

  double precision function nuclearStarClusterCrossingTimescale(radius, massStellar) result(crossingTimescale)
    !!{
    Computes the crossing timescale for a nuclear star cluster. This do not account for gas, the velocity is 
    computed internally here. 
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in   ) :: radius 
    double precision, intent(in   ) :: massStellar
    double precision                :: velocity

    if (radius<=0.0d0) then
      crossingTimescale=0.0d0
    else
      velocity         =+sqrt(                                &
        &                     +gravitationalConstant_internal &
        &                     *massStellar                    &
        &                     /radius                         &
        &                    )
      crossingTimescale=+radius/velocity
    end if
    return
  end function nuclearStarClusterCrossingTimescale

  double precision function nuclearStarClusterRelaxationTimescale(radius,massStellar,massGas,includeGasPotential) result(relaxationTime)
    !!{
      Computes the relaxation timescale for a nuclear star cluster. Here it is possible to include the gas potential in the calculations.
    !!}
    use :: Error, only : Error_Report
    implicit none 
    double precision, intent(in   )           :: radius
    double precision, intent(in   )           :: massStellar
    double precision, intent(in   ), optional :: massGas
    logical         , intent(in   ), optional :: includeGasPotential
    double precision, parameter               :: gamma=0.4d0
    double precision                          :: q                  , numberOfStars, &
      &                                          crossingTimescale
    
    ! Validate parameters.
    if (includeGasPotential.and..not.present(massGas)) &
        &   call Error_Report('The gas mass of the nuclear star cluster must be provided'//{introspection:location})
    
    if (massStellar > 0.0d0) then
      ! Let's get the number of stars from the Kroupa IMF.
      numberOfStars     = nuclearStarClusterKroupaNumberOfStars(        massStellar)
      ! Compute the crossing timescale.
      crossingTimescale = nuclearStarClusterCrossingTimescale  (radius, massStellar)
      ! Compute the relaxation timescale.
      relaxationTime =+           0.138d0 &
       &              *     numberOfStars &
       &              /log(         gamma &
       &                   *numberOfStars &
       &                  )               &
       &              *crossingTimescale
      if (includeGasPotential) then
        ! We need to correct the previous value computed to take into account the gas potential.
        ! Let's compute the factor q.
        q = massGas/massStellar
        ! Correct the relaxation time
        relaxationTime = +relaxationTime &
          &              *  (1+q)**3.0d0
      end if
    end if 
    return 
  end function nuclearStarClusterRelaxationTimescale

  double precision function nuclearStarClusterCoreCollapseTimescale(radius,massStellar,massGas,includeGasPotential) result(coreCollapseTimescale)
    !!{
      Computes the relaxation timescale for a nuclear star cluster.
    !!}
    implicit none 
    double precision, intent(in   )           :: radius
    double precision, intent(in   )           :: massStellar
    double precision, intent(in   ), optional :: massGas
    logical         , intent(in   ), optional :: includeGasPotential
   
    coreCollapseTimescale = 0.2d0*nuclearStarClusterRelaxationTimescale(radius,massStellar,massGas,includeGasPotential)
    return 
  end function nuclearStarClusterCoreCollapseTimescale

end module Nuclear_Star_Clusters_Utilities