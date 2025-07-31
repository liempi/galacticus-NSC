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
   <description>A galactic filter class which is the ``blackHoleMulti'' combination of a set of other black hole seeds.</description>
   <linkedList type="seedsList" variable="blackHoleSeeds" next="next" object="blackHoleSeeds_" objectType="blackHoleSeedsClass"/>
  </blackHoleSeeds>
  !!]


  type, extends(blackHoleSeedsClass) :: blackHoleSeedsMulti
     !!{
     A galactic filter class which is the ``blackHoleMulti'' combination of a set of other blackHoleSeeds.
     !!}
     private
     type (seedsList              ), pointer :: blackHoleSeeds     => null()

  contains
     final     ::                     blackHoleMultiDestructor
     procedure :: mass             => blackHoleMultiMass
     procedure :: spin             => blackHoleMultiSpin
     procedure :: redshift         => blackHoleMultiRedshift
     procedure :: formationChannel => blackHoleMultiFormationChannel
  end type blackHoleSeedsMulti

  interface blackHoleSeedsMulti
     !!{
     Constructors for the blackHoleMulti galactic filter class.
     !!}
     module procedure blackHoleMultiConstructorParameters
     module procedure blackHoleMultiConstructorInternal
  end interface blackHoleSeedsMulti

  ! Module-scope variable used in root finding.
  type(enumerationBlackHoleFormationChannelType) :: channel_
  !$omp threadprivate(channel_)

contains

  function blackHoleMultiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleSeedsMulti} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (blackHoleSeedsMulti    )                :: self
    type   (inputParameters        ), intent(inout) :: parameters
    type   (seedsList              ), pointer       :: blackHoleSeeds_
    integer                                         :: i

    self   %blackHoleSeeds => null()
    blackHoleSeeds_        => null()
    do i=1,parameters%copiesCount('blackHoleSeeds',zeroIfNotPresent=.true.)
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
    Internal constructor for the \refClass{blackHoleSeedsMulti} filter class.
    !!}
    implicit none
    type (blackHoleSeedsMulti    )                         :: self
    type (seedsList              ), intent(in   ), target  :: blackHoleSeeds
    type (seedsList              ),                pointer :: blackHoleSeeds_

    self       %blackHoleSeeds => blackHoleSeeds
    blackHoleSeeds_             => blackHoleSeeds
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
    Destructor for the \refClass{blackHoleSeedsMulti} galactic filter class.
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

  double precision function blackHoleMultiMass(self,node) result(mass)
   !!{
      Compute the nuclear star cluster collapse condition.
   !!}
    use :: Galacticus_Nodes, only : treeNode  
    implicit none
    class  (blackHoleSeedsMulti), intent(inout) :: self
    type   (treeNode           ), intent(inout) :: node
    type   (seedsList          ), pointer       :: blackHoleSeeds_
    logical                                     :: blackHoleMultiPasses
    ! ML: Here we need to define the logic to decide which channel is more important.
    ! In principle, now it works as first come, first served. 

    ! AB: Maybe it makes sense for each blackHoleSeeds class to also return a timescale
    ! for seed formation, and then keep whichever has the shorter timescale?
    ! Or perhaps you just sum the masses over all classes and use that as a the seed mass
    mass                 =0.0d0
    blackHoleMultiPasses =  .false.
    blackHoleSeeds_  => self%blackHoleSeeds
    do while (associated(blackHoleSeeds_))
       mass = blackHoleSeeds_%blackHoleSeeds_%mass(node)
       if (mass > 0.0d0) then 
          blackHoleMultiPasses=.true.
       end if
       if (blackHoleMultiPasses) then
          blackHoleSeeds_ => null()
       else
          blackHoleSeeds_ => blackHoleSeeds_%next
       end if
    end do
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

  double precision function blackHoleMultiRedshift(self,node) result(redshift)
    !!{
    Compute the formation redshift of the seed black hole.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode  
    implicit none
    class           (blackHoleSeedsMulti), intent(inout) :: self
    type            (treeNode           ), intent(inout) :: node
    type            (seedsList          ), pointer       :: blackHoleSeeds_
    double precision                                     :: mass
    logical                                              :: blackHoleMultiPasses
    ! ML: Here we need to define the logic to decide which channel is more important.
    ! In principle, now it works as first come, first served. 

    ! AB: Maybe it makes sense for each blackHoleSeeds class to also return a timescale
    ! for seed formation, and then keep whichever has the shorter timescale?
    ! Or perhaps you just sum the masses over all classes and use that as a the seed mass
    mass                 =0.0d0
    redshift             = 0.0d0
    blackHoleMultiPasses =  .false.
    blackHoleSeeds_  => self%blackHoleSeeds
    do while (associated(blackHoleSeeds_))
       mass = blackHoleSeeds_%blackHoleSeeds_%mass(node)
       if (mass > 0.0d0) then 
          blackHoleMultiPasses=.true.
          redshift            = blackHoleSeeds_%blackHoleSeeds_%redshift(node)
       end if
       if (blackHoleMultiPasses) then
          blackHoleSeeds_ => null()
       else
          blackHoleSeeds_ => blackHoleSeeds_%next
       end if
    end do  
   end function blackHoleMultiRedshift

  function blackHoleMultiFormationChannel (self,node) result(channel)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    type            (enumerationBlackHoleFormationChannelType)                :: channel
    class           (blackHoleSeedsMulti                     ), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    type            (seedsList                               ), pointer       :: blackHoleSeeds_
    double precision                                                          :: mass
    logical                                                                   :: blackHoleMultiPasses
   
    ! Initialize with undetermined formation
    channel              =blackHoleFormationChannelUndetermined
    blackHoleMultiPasses =  .false.
    blackHoleSeeds_      => self%blackHoleSeeds
    do while (associated(blackHoleSeeds_))
       mass = blackHoleSeeds_%blackHoleSeeds_%mass(node)
       if (mass > 0.0d0) then 
          blackHoleMultiPasses=.true.
          channel=blackHoleSeeds_%blackHoleSeeds_%formationChannel(node)
       end if
       if (blackHoleMultiPasses) then
          blackHoleSeeds_ => null()
       else
          blackHoleSeeds_ => blackHoleSeeds_%next
       end if
    end do
    return
  end function blackHoleMultiFormationChannel
