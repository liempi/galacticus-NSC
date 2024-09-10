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
  Implementation of a timescale for star formation which scales with the component Crossing time.
  !!}


  !![
  <NSCTimescale name="NSCTimescaleCrossingTime">
   <description>
    A timescale class in which the crossing timescale is computed as
    time. Specifically:
    \begin{equation}
     t_{\rm cross} = \epsilon_r \left( {R \over V(\epsilon_r)},
    \end{equation}
    where $\epsilon_r=${\normalfont \ttfamily [efficiency]}  and $r$ and $V$
    are the characteristic radius and velocity respectively of the component.
   </description>
  </NSCTimescale>
  !!]
  type, extends(NSCTimescaleClass) :: NSCTimescaleCrossingTime
     !!{
     Implementation of a timescale for star formation which scales with the Crossing time.
     !!}
     private
     double precision                                 :: efficiency      
   contains
     procedure :: timescale => CrossingTimeTimescale
  end type NSCTimescaleCrossingTime

  interface NSCTimescaleCrossingTime
     !!{
     Constructors for the {\normalfont \ttfamily CrossingTime} timescale for star formation class.
     !!}
     module procedure CrossingTimeConstructorParameters
     module procedure CrossingTimeConstructorInternal
  end interface NSCTimescaleCrossingTime

contains

  function CrossingTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily CrossingTime} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (NSCTimescaleCrossingTime)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: efficiency  

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The efficiency of star formation for the Crossing time method.</description>
      <source>parameters</source>
    </inputParameter>

    !!]
    self=NSCTimescaleCrossingTime(efficiency)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function CrossingTimeConstructorParameters

  function CrossingTimeConstructorInternal(efficiency) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily CrossingTime} timescale for star formation class.
    !!}
    implicit none
    type            (NSCTimescaleCrossingTime)                        :: self
    double precision                           , intent(in   )         :: efficiency 
    !![
    <constructorAssign variables="efficiency"/>
    !!]

    return
  end function CrossingTimeConstructorInternal

  double precision function CrossingTimeTimescale(self,node)
    !!{
    Returns the crossing timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component}. The timescale is given by

    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, treeNode
    implicit none
    class           (NSCTimescaleCrossingTime), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: treeNode
    class           (nodeComponentNSC        ), pointer       :: NSC  
    double precision                                          :: velocity     , radius
    ! Check for zero velocity.



    NSC      => node%                  NSC  ()
    radius   =  self%efficiency*NSC%radius  () !Mpc
    velocity =                  NSC%velocity()                                                                             

    if (velocity <= 0.0d0) then
       ! No well defined answer in this case.
       CrossingTimeTimescale=0.0d0
    else if (self%efficiency == 0.0d0) then
       ! No star formation occurs if the efficiency is zero.
       CrossingTimeTimescale=0.0d0
    else
       ! Get the Crossing time in Gyr.
       CrossingTimeTimescale=+Mpc_per_km_per_s_To_Gyr &
            &        *radius                  &
            &        /velocity
    end if
    return
  end function CrossingTimeTimescale
