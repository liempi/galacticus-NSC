!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) blackHoleMulti later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT blackHoleMulti WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{
Implements a galactic filter class which is the ``blackHoleMulti'' combination of a set of other blackHoleSeeds.
!!}

  !![
  <blackHoleSeeds name="blackHoleSeedsMulti">
   <description>A multi black hole seed class which is the ``blackHoleMulti'' combination of a set of other black hole seed classes.</description>
   <linkedList type="seedsList" variable="blackHoleSeeds" next="next" object="blackHoleSeeds_" objectType="blackHoleSeedsClass"/>
  </blackHoleSeeds>
  !!]

  type, extends(blackHoleSeedsClass) :: blackHoleSeedsMulti
     !!{
     A black hole seeds class which is the ``blackHoleMulti'' combination of a set of other blackHoleSeeds.
     !!}
     private
     type   (seedsList), pointer :: blackHoleSeeds              => null()
     ! This variable is used in order to determinate the size of the array that will store the results of the
     ! blackHoleMultiMass result for each class. We should have some sort of method that detects if the is more
     ! than one class that returns a non zero mass for the black hole seed at the same timestep.
     ! The idea is to create a new method that will check that, if that is true, then we should obtain the index
     ! of the element that has the shorter timescale and use as the prefered seeding mechanism.
     integer                     :: blackHoleSeedsClassCounts_

  contains
     !![
     <methods>
       <method description="Returns the index of the shortest timescale in order to return the mass of the black hole seed." method="index"/>
     </methods>
     !!]
     final     ::                     blackHoleMultiDestructor
     procedure :: index            => blackHoleMultiIndex 
     procedure :: timescale        => blackHoleMultiTimescale
     procedure :: mass             => blackHoleMultiMass
     procedure :: spin             => blackHoleMultiSpin
     procedure :: formationChannel => blackHoleMultiFormationChannel
  end type blackHoleSeedsMulti

  interface blackHoleSeedsMulti
     !!{
     Constructors for the blackHoleMulti seed class.
     !!}
     module procedure blackHoleMultiConstructorParameters
     module procedure blackHoleMultiConstructorInternal
  end interface blackHoleSeedsMulti

  ! Module-scope variable used in root finding.
  type            (enumerationBlackHoleFormationChannelType), dimension(:), allocatable :: channel_
  double precision                                          , dimension(:), allocatable :: masses_  , timescales_
  !$omp threadprivate(channel_, masses_, timescales_)

contains

  function blackHoleMultiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleSeedsMulti} seed class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (blackHoleSeedsMulti)                :: self
    type   (inputParameters    ), intent(inout) :: parameters
    type   (seedsList          ), pointer       :: blackHoleSeeds_
    integer                                     :: i

    self%blackHoleSeeds            => null()
    blackHoleSeeds_                => null()
    self%blackHoleSeedsClassCounts_=  parameters%copiesCount('blackHoleSeeds',zeroIfNotPresent=.true.)

    do i=1,self%blackHoleSeedsClassCounts_
       if (associated(blackHoleSeeds_)) then
          allocate(blackHoleSeeds_%next)
          blackHoleSeeds_ => blackHoleSeeds_%next
       else
          allocate(self%blackHoleSeeds)
          blackHoleSeeds_ => self%blackHoleSeeds
       end if
      !![
       <objectBuilder class="blackHoleSeeds" name="blackHoleSeeds_%blackHoleSeeds_" source="parameters" copy="i" />
      !!]
    end do

    !![
    <inputParametersValidate source="parameters" multiParameters="blackHoleSeeds"/>
    !!]
    return
  end function blackHoleMultiConstructorParameters

  function blackHoleMultiConstructorInternal(blackHoleSeeds) result(self)
    !!{
    Internal constructor for the \refClass{blackHoleSeedsMulti} seed class.
    !!}
    implicit none
    type (blackHoleSeedsMulti)                         :: self
    type (seedsList          ), intent(in   ), target  :: blackHoleSeeds
    type (seedsList          ),                pointer :: blackHoleSeeds_

    self           %blackHoleSeeds => blackHoleSeeds
    blackHoleSeeds_                => blackHoleSeeds
    do while (associated(blackHoleSeeds_))
       !![
       <referenceCountIncrement owner="blackHoleSeeds_" object="blackHoleSeeds_"/>
       !!]
       blackHoleSeeds_ => blackHoleSeeds_%next
    end do
    return
  end function blackHoleMultiConstructorInternal

  subroutine blackHoleMultiDestructor(self)
    !!{
    Destructor for the \refClass{blackHoleSeedsMulti} seed class.
    !!}
    implicit none
    type(blackHoleSeedsMulti), intent(inout) :: self
    type(seedsList          ), pointer       :: blackHoleSeeds_, blackHoleSeedsNext

    if (associated(self%blackHoleSeeds)) then
       blackHoleSeeds_ => self%blackHoleSeeds
       do while (associated(blackHoleSeeds_))
          blackHoleSeedsNext => blackHoleSeeds_%next
          !![
          <objectDestructor name="blackHoleSeeds_%blackHoleSeeds_"/>
          !!]
          deallocate(blackHoleSeeds_)
          blackHoleSeeds_ => blackHoleSeedsNext
       end do
    end if
    return
  end subroutine blackHoleMultiDestructor

  integer function blackHoleMultiIndex(self,node)
   !!{
   Return the index of the shortest timescale.
   !!}
   use :: Galacticus_Nodes, only : treeNode
   implicit none
   class(blackHoleSeedsMulti), intent(inout) :: self
   type (treeNode           ), intent(inout) :: node
   !$GLC attributes unused :: self, node

   blackHoleMultiIndex=minloc(timescales_,1)
   return
  end function blackHoleMultiIndex

  double precision function blackHoleMultiTimescale(self,node)
   !!{
      Return the black hole seed timescale according to condition.
   !!}
    use :: Galacticus_Nodes, only : treeNode  
    implicit none
    class  (blackHoleSeedsMulti), intent(inout) :: self
    type   (treeNode           ), intent(inout) :: node

    blackHoleMultiTimescale=timescales_(self%index(node))
    return
  end function blackHoleMultiTimescale

  double precision function blackHoleMultiMass(self,node) result(mass)
   !!{
      Compute the black hole masses according to condition.
   !!}
    use :: Galacticus_Nodes, only : treeNode  
    implicit none
    class           (blackHoleSeedsMulti), intent(inout) :: self
    type            (treeNode           ), intent(inout) :: node
    type            (seedsList          ), pointer       :: blackHoleSeeds_
    double precision                                     :: timescale
    integer                                              :: i               , minimumTimescaleIndex

    i                =0
    mass             =0.0d0
    timescale        =0.0d0
    blackHoleSeeds_  => self%blackHoleSeeds
    do while (associated(blackHoleSeeds_))
       masses_    (i)= blackHoleSeeds_%blackHoleSeeds_%            mass(node)
       timescales_(i)= blackHoleSeeds_%blackHoleSeeds_%       timescale(node)
       channel_   (i)= blackHoleSeeds_%blackHoleSeeds_%formationChannel(node)
       i = i+1
    end do

    mass = masses_(minimumTimescaleIndex)

    return
  end function blackHoleMultiMass

  double precision function blackHoleMultiSpin(self,node) result(spin)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsMulti), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Assume zero spin.
    spin=0.0d0
    return
  end function blackHoleMultiSpin

  function blackHoleMultiFormationChannel (self,node) result(channel)
    !!{
    Retuns the channel formation of the seed black hole.
    !!}
    implicit none
    type            (enumerationBlackHoleFormationChannelType)                :: channel
    class           (blackHoleSeedsMulti                     ), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    channel=channel_(self%index(node))
    return
  end function blackHoleMultiFormationChannel
