!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which interfaces with the Local Group database.

module Interface_Local_Group_DB
  !% Interfaces with the Local Group database.
  use FoX_DOM, only : node, nodeList
  private
  public  :: localGroupDB, vector3D

  !# <generic identifier="Type">
  !#  <instance label="Double" intrinsic="double precision"                />
  !#  <instance label="VarStr" intrinsic="type            (varying_string)"/>
  !#  <instance label="Vector" intrinsic="type            (vector3D      )"/>
  !# </generic>
  
  !# <enumeration>
  !#  <name>comparison</name>
  !#  <description>Comparison operators.</description>
  !#  <visibility>public</visibility>
  !#  <entry label="equals"     />
  !#  <entry label="greaterThan"/>
  !#  <entry label="lessThan"   />
  !# </enumeration>
  
  !# <enumeration>
  !#  <name>setOperator</name>
  !#  <description>Set operators.</description>
  !#  <visibility>public</visibility>
  !#  <entry label="intersection"/>
  !#  <entry label="union"/>
  !#  <entry label="relativeComplement"/>
  !# </enumeration>

  type :: vector3D
     !% Vector type.
     double precision, dimension(3) :: x
   contains
     procedure ::                 vector3DEquals
     procedure ::                 vector3DComparisonUnimplemented
     generic   :: operator(==) => vector3DEquals
     generic   :: operator(<)  => vector3DComparisonUnimplemented
     generic   :: operator(>)  => vector3DComparisonUnimplemented
  end type vector3D
  
  type :: localGroupDB
     !% Local Group database class.
     private
     type   (node    ), pointer                   :: database
     type   (nodeList), pointer                   :: galaxies
     logical          , allocatable, dimension(:) :: selected
   contains
     !@ <objectMethods>
     !@   <object>localGroupDB</object>
     !@   <objectMethod>
     !@     <method>getProperty</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater}name\argin, \textcolor{red}{\textless type(varying\_string)(:)\textgreater}property\arginout</arguments>
     !@     <description>Return an array of values of the named property for the current selection.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>select</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater}name\argin, \textcolor{red}{\textless character(len=*)\textgreater}value\argin</arguments>
     !@     <description>Select all galaxies in the current selection where the named property has the given value.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>selectAll</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Select all galaxies in the database.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>update</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Update the database.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                localGroupDBDestructor
     procedure ::                localGroupDBGetProperty{Type¦label}
     generic   :: getProperty => localGroupDBGetProperty{Type¦label}
     procedure ::                localGroupDBSelect{Type¦label}
     generic   :: select      => localGroupDBSelect{Type¦label}
     procedure :: selectAll   => localGroupDBSelectAll
     procedure :: update      => localGroupDBUpdate
  end type localGroupDB

  interface localGroupDB
     !% Constructors for the Local Group database class.
     module procedure localGroupDBConstructorInternal
  end interface localGroupDB
  
contains

  function localGroupDBConstructorInternal() result(self)
    !% Constructor for the Local Group database class.
    use FoX_DOM           , only : getElementsByTagName, getLength
    use IO_XML            , only : XML_Parse           , XML_Get_First_Element_By_Tag_Name
    use Galacticus_Paths  , only : galacticusPath      , pathTypeDataStatic
    use ISO_Varying_String, only : char
    implicit none
    type(localGroupDB) :: self

    self%database => XML_Parse(char(galacticusPath(pathTypeDataStatic))//"observations/localGroup/localGroupSatellites.xml")
    self%galaxies => getElementsByTagName(XML_Get_First_Element_By_Tag_Name(self%database,'galaxies'),'galaxy')
    allocate(self%selected(0:getLength(self%galaxies)-1))
    self%selected=.false.
    return
  end function localGroupDBConstructorInternal
  
  subroutine localGroupDBDestructor(self)
    !% Destructor for the Local Group database class.
    use FoX_DOM, only : destroy
    implicit none
    type(localGroupDB), intent(inout) :: self

    call destroy(self%database)
    return
  end subroutine localGroupDBDestructor
  
  subroutine localGroupDBGetProperty{Type¦label}(self,name,property,isPresent)
    !% Get a named text property from the Local Group database.
    use FoX_DOM           , only : getLength              , item          , hasAttribute      , getElementsByTagName, &
         &                         getAttributeNode       , getTextContent, extractDataContent
    use Galacticus_Error  , only : Galacticus_Error_Report
    use ISO_Varying_String
    implicit none
    class           (localGroupDB  ), intent(inout)                                      :: self
    character       (len=*         ), intent(in   )                                      :: name
    {Type¦intrinsic}                , intent(inout), allocatable, dimension(:)           :: property
    logical                         , intent(inout), allocatable, dimension(:), optional :: isPresent
    type            (node          )               , pointer                             :: galaxy    , propertyElement, &
         &                                                                                  attribute
    type            (nodeList      )               , pointer                             :: properties
    integer                                                                              :: i         , countGalaxies
    
    ! Iterate over galaxies to count how many have this property.
    countGalaxies=0
    do i=0,getLength(self%galaxies)-1
       if (.not.self%selected(i)) cycle
       galaxy => item(self%galaxies,i)
       if (present(isPresent)) then
          countGalaxies=countGalaxies+1
       else if (hasAttribute(galaxy,name)) then
          countGalaxies=countGalaxies+1
       else
          properties => getElementsByTagName(galaxy,name)
          if      (getLength(properties) > 1) then
             call Galacticus_Error_Report('galaxy has multiple entries for named property'//{introspection:location})
          else if (getLength(properties) == 1) then
             countGalaxies=countGalaxies+1
          end if
       end if
    end do
    allocate(property(countGalaxies))
    if (present(isPresent)) then
       allocate(isPresent(countGalaxies))
       isPresent=.true.
    end if
    ! Popuate the array.
    countGalaxies=0
    do i=0,getLength(self%galaxies)-1
       if (.not.self%selected(i)) cycle
       galaxy => item(self%galaxies,i)
       if (hasAttribute(galaxy,name)) then
          countGalaxies   =  countGalaxies+1
          attribute => getAttributeNode(galaxy,name)
       else
          properties => getElementsByTagName(galaxy,name)
          if (getLength(properties) == 1) then
             countGalaxies=countGalaxies+1
             propertyElement => item(properties,0)
             attribute       => getAttributeNode(propertyElement,'value')
          else
             if (present(isPresent)) then
                countGalaxies=countGalaxies+1
                isPresent(countGalaxies)=.false.
                {Type¦match¦^Double$¦property(countGalaxies)=0.0d0¦}
                {Type¦match¦^Vector¦property(countGalaxies)%x=0.0d0¦}
                {Type¦match¦^VarStr$¦property(countGalaxies)='not present'¦}
             else
                attribute       => getAttributeNode(galaxy,'name')
                call Galacticus_Error_Report('property "'//name//'" is not present in galaxy "'//getTextContent(attribute)//'"'//{introspection:location})
             end if
          end if
       end if
       {Type¦match¦^Double$¦call extractDataContent(attribute,property(countGalaxies))¦}
       {Type¦match¦^Vector¦call extractDataContent(attribute,property(countGalaxies)%x)¦}
       {Type¦match¦^VarStr$¦property(countGalaxies)=getTextContent(attribute)¦}
    end do
    return
  end subroutine localGroupDBGetProperty{Type¦label}
    
  subroutine localGroupDBSelectAll(self)
    !% Select all galaxies in the database.
    implicit none
    class(localGroupDB), intent(inout) :: self
    
    self%selected=.true.
    return
  end subroutine localGroupDBSelectAll

  subroutine localGroupDBSelect{Type¦label}(self,name,value,comparison,setOperator)
    !% Impose a selection on the database.
    use ISO_Varying_String
    use FoX_DOM           , only : getLength              , getAttributeNode, getTextContent, item
    use Galacticus_Error  , only : Galacticus_Error_Report
    implicit none
    class           (localGroupDB  ), intent(inout)               :: self
    character       (len=*         ), intent(in   )               :: name
    {Type¦intrinsic}                , intent(in)                  :: value
    integer                         , intent(in   )               :: comparison      , setOperator
    {Type¦intrinsic}                , allocatable  , dimension(:) :: values
    logical                         , allocatable  , dimension(:) :: selectedCurrent , isPresent
    type            (node          )               , pointer      :: galaxy          , attribute
    integer                                                       :: i
    logical                                                       :: comparisonResult

    allocate(selectedCurrent(0:getLength(self%galaxies)-1))
    selectedCurrent=self%selected
    self%selected=.true.
    call self%getProperty(name,values,isPresent)
    do i=0,getLength(self%galaxies)-1
       if (.not.selectedCurrent(i) .and. (setOperator == setOperatorIntersection .or.  setOperator == setOperatorRelativeComplement)) cycle
       if (.not.isPresent(i+1)) then
          galaxy => item(self%galaxies,i)
          attribute       => getAttributeNode(galaxy,'name')
         call Galacticus_Error_Report('property "'//name//'" is not present in selected galaxy "'//getTextContent(attribute)//'"'//{introspection:location})
       end if
       select case (comparison)
       case (comparisonEquals     )
          comparisonResult=values(i+1) == value
       case (comparisonLessThan   )
          comparisonResult=values(i+1) <  value
       case (comparisonGreaterThan)
          comparisonResult=values(i+1) >  value
       case default
          comparisonResult=.false.
          call Galacticus_Error_Report('unknown comparison operator'//{introspection:location})
       end select       
       select case (setOperator)
       case (setOperatorIntersection)
          selectedCurrent(i)=selectedCurrent(i) .and. comparisonResult
       case (setOperatorUnion)
          selectedCurrent(i)=selectedCurrent(i) .or. comparisonResult
       case (setOperatorRelativeComplement)
          selectedCurrent(i)=selectedCurrent(i) .and. .not.comparisonResult
       case default
          call Galacticus_Error_Report('unknown set operator'       //{introspection:location})
       end select
    end do
    self%selected=selectedCurrent
    return
  end subroutine localGroupDBSelect{Type¦label}

  subroutine localGroupDBUpdate(self)
    !% Update the database.
    use FoX_DOM                         , only : getLength          , item               , getElementsByTagName, getAttributeNode, &
         &                                       extractDataContent , createElementNS    , getNamespaceURI     , appendChild     , &
         &                                       setAttribute       , serialize          , hasAttribute        , getTextContent
    use Numerical_Constants_Astronomical, only : arcminutesToDegrees, arcsecondsToDegrees, hoursToDegrees      , minutesToDegrees, &
         &                                       secondsToDegrees   , degreesToRadians
    use Galacticus_Paths                , only : galacticusPath     , pathTypeDataStatic
    use ISO_Varying_String              , only : char
    implicit none
    class           (localGroupDB  ), intent(inout) :: self
    type            (node          ), pointer       :: galaxy                        , attribute                       , newNode
    type            (nodeList      ), pointer       :: propertyList1                 , propertyList2                   , propertyList3                   , propertyList4
    character       (len=32        ), dimension(4)  :: uncertainties
    double precision                , dimension(3)  :: positionMilkyWay              , positionM31                     , position                        , uncertaintyPosition  , &
         &                                             uncertaintyPositionMilkyWay   , uncertaintyPositionM31
    integer                                         :: i                             , j                               , indexMilkyWay                   , indexM31
    logical                                         :: hasSexagesimal                , hasDecimal                      , hasHeliocentric                 , hasModulus           , &
         &                                             hasDistance                   , hasDeclination                  , hasRightAscension               , hasPosition
    double precision                                :: declinationSexagesimalDegrees , declinationSexagesimalArcMinutes, declinationSexagesimalArcSeconds, declinationDecimal   , &
         &                                             rightAscensionSexagesimalHours, rightAscensionSexagesimalMinutes, rightAscensionSexagesimalSeconds, rightAscensionDecimal, &
         &                                             distanceHeliocentric          , distanceModulus                 , uncertainty1                    , uncertainty2         , &
         &                                             positionHeliocentricX         , positionHeliocentricY           , positionHeliocentricZ           , distance
    character       (len=64        )                :: textContent

    ! Set names of uncertainties.
    uncertainties(1)='uncertainty'
    uncertainties(2)='uncertaintyLow'
    uncertainties(3)='uncertaintyHigh'
    uncertainties(4)='uncertaintySystematic'
    ! Identify Milky Way and M31.
    do i=0,getLength(self%galaxies)-1
       galaxy    => item            (self%galaxies,i     )
       attribute => getAttributeNode(     galaxy  ,'name')
       if (getTextContent(attribute) == "The Galaxy") indexMilkyWay=i
       if (getTextContent(attribute) == "Andromeda" ) indexM31     =i
    end do
    ! Iterate over galaxies.
    do i=0,getLength(self%galaxies)-1
       galaxy => item(self%galaxies,i)
       ! Look for declinations.
       propertyList1  => getElementsByTagName(galaxy,'declinationSexagesimal')
       hasSexagesimal =  getLength(propertyList1) == 1
       propertyList2  => getElementsByTagName(galaxy,'declinationDecimal'    )
       hasDecimal     =  getLength(propertyList2) == 1
       if (hasSexagesimal .and. .not. hasDecimal) then
          ! Convert sexagesimal declination to decimal.
          attribute => getAttributeNode(item(propertyList1,0),'degrees'  )
          call extractDataContent(attribute,declinationSexagesimalDegrees   )
          attribute => getAttributeNode(item(propertyList1,0),'arcminutes')
          call extractDataContent(attribute,declinationSexagesimalArcMinutes)
          attribute => getAttributeNode(item(propertyList1,0),'arcseconds')
          call extractDataContent(attribute,declinationSexagesimalArcSeconds)
          declinationDecimal=+     declinationSexagesimalDegrees                                                       &
               &             +sign(declinationSexagesimalArcMinutes,declinationSexagesimalDegrees)*arcminutesToDegrees &
               &             +sign(declinationSexagesimalArcSeconds,declinationSexagesimalDegrees)*arcsecondsToDegrees
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'declinationDecimal')
          write (textContent,'(f16.12)') declinationDecimal
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       else if (.not.hasSexagesimal .and. hasDecimal) then
          ! Convert decimal declination to sexagesimal.
          attribute => getAttributeNode(item(propertyList2,0),'value')
          call extractDataContent(attribute,declinationDecimal)
          declinationSexagesimalDegrees   =int( declinationDecimal                                                   )
          declinationSexagesimalArcminutes=int((declinationDecimal-declinationSexagesimalDegrees)/arcminutesToDegrees)
          declinationSexagesimalArcseconds=    (declinationDecimal-declinationSexagesimalDegrees-declinationSexagesimalArcminutes*arcminutesToDegrees)/arcsecondsToDegrees
          declinationSexagesimalArcminutes=abs(declinationSexagesimalArcminutes)
          declinationSexagesimalArcseconds=abs(declinationSexagesimalArcseconds)
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'declinationSexagesimal')
          write (textContent,'(f16.12)') declinationSexagesimalDegrees
          call setAttribute(newNode,"degrees"   ,trim(adjustl(textContent)))
          write (textContent,'(f16.12)') declinationSexagesimalArcMinutes
          call setAttribute(newNode,"arcminutes",trim(adjustl(textContent)))
          write (textContent,'(f16.12)') declinationSexagesimalArcSeconds
          call setAttribute(newNode,"arcseconds",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       end if
       ! Look for right ascensions.
       propertyList1 => getElementsByTagName(galaxy,'rightAscensionSexagesimal')
       hasSexagesimal=getLength(propertyList1) == 1
       propertyList2 => getElementsByTagName(galaxy,'rightAscensionDecimal'    )
       hasDecimal=getLength(propertyList2) == 1
       if (hasSexagesimal .and. .not. hasDecimal) then
          ! Convert sexagesimal right ascension to decimal.
          attribute => getAttributeNode(item(propertyList1,0),'hours'  )
          call extractDataContent(attribute,rightAscensionSexagesimalHours  )
          attribute => getAttributeNode(item(propertyList1,0),'minutes')
          call extractDataContent(attribute,rightAscensionSexagesimalMinutes)
          attribute => getAttributeNode(item(propertyList1,0),'seconds')
          call extractDataContent(attribute,rightAscensionSexagesimalSeconds)
          rightAscensionDecimal=+rightAscensionSexagesimalHours  *hoursToDegrees   &
               &                +rightAscensionSexagesimalMinutes*minutesToDegrees &
               &                +rightAscensionSexagesimalSeconds*secondsToDegrees
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'rightAscensionDecimal')
          write (textContent,'(f16.12)') rightAscensionDecimal
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       else if (.not.hasSexagesimal .and. hasDecimal) then
          ! Convert decimal right ascension to sexagesimal.
          attribute => getAttributeNode(item(propertyList2,0),'value')
          call extractDataContent(attribute,rightAscensionDecimal)
          rightAscensionSexagesimalHours  =int( rightAscensionDecimal/hoursToDegrees                                                 )
          rightAscensionSexagesimalMinutes=int((rightAscensionDecimal-rightAscensionSexagesimalHours*hoursToDegrees)/minutesToDegrees)
          rightAscensionSexagesimalSeconds=(rightAscensionDecimal-rightAscensionSexagesimalHours*hoursToDegrees-rightAscensionSexagesimalMinutes*minutesToDegrees)/secondsToDegrees
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'rightAscensionSexagesimal')
          write (textContent,'(f16.12)') rightAscensionSexagesimalHours
          call setAttribute(newNode,"hours"  ,trim(adjustl(textContent)))
          write (textContent,'(f16.12)') rightAscensionSexagesimalMinutes
          call setAttribute(newNode,"minutes",trim(adjustl(textContent)))
          write (textContent,'(f16.12)') rightAscensionSexagesimalSeconds
          call setAttribute(newNode,"seconds",trim(adjustl(textContent)))
          newNode => appendChild(galaxy,newNode)
       end if
       ! Look for distances
       propertyList1 => getElementsByTagName(galaxy,'distanceHeliocentric')
       hasHeliocentric=getLength(propertyList1) == 1
       propertyList2 => getElementsByTagName(galaxy,'distanceModulus'     )
       hasModulus=getLength(propertyList2) == 1
       if (hasHeliocentric .and. .not. hasModulus) then
          ! Convert heliocentric distance to a distance modulus.
          attribute => getAttributeNode(item(propertyList1,0),'value')
          call extractDataContent(attribute,distanceHeliocentric)
          distanceModulus=25.0d0+5.0d0*log10(distanceHeliocentric)
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'distanceModulus')
          write (textContent,'(f16.12)') distanceModulus
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(item(propertyList1,0),trim(uncertainties(j)))) then
                attribute => getAttributeNode(item(propertyList1,0),trim(uncertainties(j)))
                call extractDataContent(attribute,uncertainty1)
                uncertainty2=5.0d0/log(10.0d0)*uncertainty1/distanceHeliocentric
                write (textContent,'(f16.12)') uncertainty2
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       else if (.not.hasHeliocentric .and. hasModulus) then
          ! Convert distance modules to a heliocentric distance.          
          attribute => getAttributeNode(item(propertyList2,0),'value')
          call extractDataContent(attribute,distanceModulus)
          distanceHeliocentric=10.0d0**((distanceModulus-25.0d0)/5.0d0)
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'distanceHeliocentric')
          write (textContent,'(f16.12)') distanceHeliocentric
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(item(propertyList2,0),trim(uncertainties(j)))) then
                attribute => getAttributeNode(item(propertyList2,0),trim(uncertainties(j)))
                call extractDataContent(attribute,uncertainty1)
                uncertainty2=distanceHeliocentric*log(10.0d0)/5.0d0
                write (textContent,'(f16.12)') uncertainty2
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       end if
       ! Compute Cartesian heliocentric coordinates.
       propertyList1 => getElementsByTagName(galaxy,'distanceHeliocentric'         )
       hasDistance=getLength(propertyList1) == 1
       propertyList2 => getElementsByTagName(galaxy,'rightAscensionDecimal'        )
       hasRightAscension=getLength(propertyList2) == 1
       propertyList3 => getElementsByTagName(galaxy,'declinationDecimal'           )
       hasDeclination=getLength(propertyList3) == 1
       propertyList4 => getElementsByTagName(galaxy,'positionHeliocentricCartesian')
       hasPosition=getLength(propertyList4) == 1
       if (hasDistance .and. hasRightAscension .and. hasDeclination .and. .not. hasPosition) then
          attribute => getAttributeNode(item(propertyList1,0),'value')
          call extractDataContent(attribute,distanceHeliocentric )         
          attribute => getAttributeNode(item(propertyList2,0),'value')
          call extractDataContent(attribute,rightAscensionDecimal)         
          attribute => getAttributeNode(item(propertyList3,0),'value')
          call extractDataContent(attribute,declinationDecimal   )
          positionHeliocentricX=distanceHeliocentric*cos(declinationDecimal*degreesToRadians)*cos(rightAscensionDecimal*degreesToRadians)
          positionHeliocentricY=distanceHeliocentric*cos(declinationDecimal*degreesToRadians)*sin(rightAscensionDecimal*degreesToRadians)
          positionHeliocentricZ=distanceHeliocentric*sin(declinationDecimal*degreesToRadians)
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'positionHeliocentricCartesian')
          write (textContent,'(e16.8,1x,e16.8,1x,e16.8)') positionHeliocentricX,positionHeliocentricY,positionHeliocentricZ
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          do j=1,size(uncertainties)
             if (hasAttribute(item(propertyList1,0),trim(uncertainties(j)))) then
                attribute => getAttributeNode(item(propertyList1,0),trim(uncertainties(j)))
                call extractDataContent(attribute,uncertainty1)
                positionHeliocentricX=abs(uncertainty1*cos(declinationDecimal*degreesToRadians)*cos(rightAscensionDecimal*degreesToRadians))
                positionHeliocentricY=abs(uncertainty1*cos(declinationDecimal*degreesToRadians)*sin(rightAscensionDecimal*degreesToRadians))
                positionHeliocentricZ=abs(uncertainty1*sin(declinationDecimal*degreesToRadians))

                write (textContent,'(e16.8,1x,e16.8,1x,e16.8)') positionHeliocentricX,positionHeliocentricY,positionHeliocentricZ
                call setAttribute(newNode,trim(uncertainties(j)),trim(adjustl(textContent)))
             end if
          end do
          newNode => appendChild(galaxy,newNode)
       end if
    end do
    ! Get heliocentric Cartesian positions of Milky Way and M31.
    galaxy        => item                (     self%galaxies   ,indexMilkyWay                  )
    propertyList1 => getElementsByTagName(          galaxy     ,'positionHeliocentricCartesian')
    attribute     => getAttributeNode    (item(propertyList1,0),'value'                        )
    call extractDataContent(attribute,positionMilkyWay           )             
    attribute => getAttributeNode(item(propertyList1,0),'uncertainty')
    call extractDataContent(attribute,uncertaintyPositionMilkyWay)             
    galaxy        => item                (     self%galaxies   ,indexM31                       )
    propertyList1 => getElementsByTagName(          galaxy     ,'positionHeliocentricCartesian')
    attribute     => getAttributeNode    (item(propertyList1,0),'value'                        )
    call extractDataContent(attribute,positionM31)             
    attribute => getAttributeNode(item(propertyList1,0),'uncertainty')
    call extractDataContent(attribute,uncertaintyPositionM31)             
    ! Compute galacto-centric distances.    
    do i=0,getLength(self%galaxies)-1
       ! Milky Way.
       galaxy        => item(self%galaxies,i)
       propertyList1 => getElementsByTagName(galaxy,'positionHeliocentricCartesian')
       hasPosition   =  getLength(propertyList1) == 1
       propertyList2 => getElementsByTagName(galaxy,'distanceMilkyWay'             )
       hasDistance   =  getLength(propertyList2) == 1
       if (hasPosition .and. .not.hasDistance) then
          attribute => getAttributeNode(item(propertyList1,0),'value')
          call extractDataContent(attribute,position)
          distance=sqrt(sum((position-positionMilkyWay)**2))
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'distanceMilkyWay')
          write (textContent,'(e16.8)') distance
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          if (hasAttribute(item(propertyList1,0),'uncertainty')) then
             attribute => getAttributeNode(item(propertyList1,0),'uncertainty')
             call extractDataContent(attribute,uncertaintyPosition)
             if (distance > 0.0d0) then
                uncertainty2=sqrt(sum((position-positionMilkyWay)**2*(uncertaintyPosition**2+uncertaintyPositionMilkyWay**2)))/distance
             else
                uncertainty2=0.0d0
             end if
             write (textContent,'(e16.8)') uncertainty2
             call setAttribute(newNode,'uncertainty',trim(adjustl(textContent)))
          end if
          newNode => appendChild(galaxy,newNode)
       end if
       ! M31.
       galaxy        => item(self%galaxies,i)
       propertyList1 => getElementsByTagName(galaxy,'positionHeliocentricCartesian')
       hasPosition   =  getLength(propertyList1) == 1
       propertyList2 => getElementsByTagName(galaxy,'distanceM31'                  )
       hasDistance   =  getLength(propertyList2) == 1
       if (hasPosition .and. .not.hasDistance) then
          attribute => getAttributeNode(item(propertyList1,0),'value')
          call extractDataContent(attribute,position)
          distance=sqrt(sum((position-positionM31)**2))
          newNode => createElementNS(self%database,getNamespaceURI(self%database),'distanceM31')
          write (textContent,'(e16.8)') distance
          call setAttribute(newNode,"value",trim(adjustl(textContent)))
          if (hasAttribute(item(propertyList1,0),'uncertainty')) then
             attribute => getAttributeNode(item(propertyList1,0),'uncertainty')
             call extractDataContent(attribute,uncertaintyPosition)
             if (distance > 0.0d0) then
                uncertainty2=sqrt(sum((position-positionM31)**2*(uncertaintyPosition**2+uncertaintyPositionM31**2)))/distance
             else
                uncertainty2=0.0d0
             end if
             write (textContent,'(e16.8)') uncertainty2
             call setAttribute(newNode,'uncertainty',trim(adjustl(textContent)))
          end if
          newNode => appendChild(galaxy,newNode)
       end if
    end do
    ! Write out the updated database.    
    call serialize(self%database,char(galacticusPath(pathTypeDataStatic))//"observations/localGroup/localGroupSatellites.xml")
    return
  end subroutine localGroupDBUpdate

  logical function vector3DEquals(self,other)
    !% Equality comparison operator for two 3D vectors.
    implicit none
    class(vector3D), intent(in   ) :: self, other

    vector3DEquals=all(self%x == other%x)
    return
  end function vector3DEquals
  
  logical function vector3DComparisonUnimplemented(self,other)
    !% Unimplemented comparison operators for 3D vectors.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(vector3D), intent(in   ) :: self, other
    !GCC$ attributes unused :: self, other
    
    vector3DComparisonUnimplemented=.false.
    call Galacticus_Error_Report('comparison operator is unimplemented'//{introspection:location})
    return
  end function vector3DComparisonUnimplemented
  
end module Interface_Local_Group_DB