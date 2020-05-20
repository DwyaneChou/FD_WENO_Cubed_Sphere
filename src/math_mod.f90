module math_mod
  use constants_mod
  implicit none
  
  integer(i4) :: OrientLeft  = 1
  integer(i4) :: OrientRight = 2
  integer(i4) :: OrientOn    = 0
  
  contains
  ! cross product for 3D vector
  function cross(a, b)
    real(r8), DIMENSION(3)             :: cross
    real(r8), DIMENSION(3), INTENT(IN) :: a, b
  
    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross
  
  function norm(a)
    real(r8)            :: norm
    real(r8),intent(in) :: a(:)
    
    norm = sqrt(dot_product(a,a))
    
  end function norm
    
  ! Conert coordinate from cartesian to sphere. The same function of which in Matlab
  subroutine cart2sph(longitude,latitude,x,y,z)
    real(r8), intent(out) :: longitude
    real(r8), intent(out) :: latitude
    real(r8), intent(in ) :: x
    real(r8), intent(in ) :: y
    real(r8), intent(in ) :: z
    
    longitude = atan2(y,x)
    latitude  = atan2(z,sqrt(x*x + y*y))
  end subroutine cart2sph
  
  ! Conert coordinate from sphere to cartesian. The same function of which in Matlab
  subroutine sph2cart(x,y,z,longitude,latitude,radius)
    real(r8), intent(out) :: x
    real(r8), intent(out) :: y
    real(r8), intent(out) :: z
    real(r8), intent(in ) :: longitude
    real(r8), intent(in ) :: latitude
    real(r8), intent(in ) :: radius
    
    x = radius * cos(longitude) * cos(latitude)
    y = radius * sin(longitude) * cos(latitude)
    z = radius * sin(latitude )
  end subroutine sph2cart
  
  ! Rotation matrix for rotations around x-axis. The same function of which in Matlab
  subroutine rotx(mtx,ang)
    real(r8),intent(out) :: mtx(3,3)
    real(r8),intent(in ) :: ang
    
    mtx(1,1) = 1._r16
    mtx(1,2) = 0._r16
    mtx(1,3) = 0._r16
    mtx(2,1) = 0._r16
    mtx(2,2) = cos(ang)
    mtx(2,3) = -sin(ang)
    mtx(3,1) = 0._r16
    mtx(3,2) = -mtx(2,3)
    mtx(3,3) = mtx(2,2)
  
  end subroutine rotx
  
  ! Rotation matrix for rotations around y-axis. The same function of which in Matlab
  subroutine roty(mtx,ang)
    real(r8),intent(out) :: mtx(3,3)
    real(r8),intent(in ) :: ang
    
    mtx(1,1) = cos(ang)
    mtx(1,2) = 0._r16
    mtx(1,3) = sin(ang)
    mtx(2,1) = 0._r16
    mtx(2,2) = 1._r16
    mtx(2,3) = 0._r16
    mtx(3,1) = -mtx(1,3)
    mtx(3,2) = 0._r16
    mtx(3,3) = mtx(1,1)
  
  end subroutine roty
  
  ! Rotation matrix for rotations around z-axis. The same function of which in Matlab
  subroutine rotz(mtx,ang)
    real(r8),intent(out) :: mtx(3,3)
    real(r8),intent(in ) :: ang
    
    mtx(1,1) = cos(ang)
    mtx(1,2) = -sin(ang)
    mtx(1,3) = 0._r16
    mtx(2,1) = -mtx(1,2)
    mtx(2,2) = mtx(1,1)
    mtx(2,3) = 0._r16
    mtx(3,1) = 0._r16
    mtx(3,2) = 0._r16
    mtx(3,3) = 1._r16
  
  end subroutine rotz
  
  subroutine wind_sph2cart(uc,vc,wc,us,vs,lon,lat)
    real,intent(out) :: uc,vc,wc
    real,intent(in ) :: us,vs
    real,intent(in ) :: lon,lat
    
    real s2c_u(3),s2c_v(3)
    
    s2c_u = [-sin(lon), cos(lon), 0.]
    s2c_v = [-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)]
    
    uc = s2c_u(1)*us + s2c_v(1)*vs
    vc = s2c_u(2)*us + s2c_v(2)*vs
    wc = s2c_u(3)*us + s2c_v(3)*vs
  
  end subroutine wind_sph2cart
  
  subroutine wind_cart2sph(us,vs,uc,vc,wc,lon,lat)
    real,intent(out) :: us,vs
    real,intent(in ) :: uc,vc,wc
    real,intent(in ) :: lon,lat
    
    real c2s_u(2),c2s_v(2),c2s_w(2)
    
    c2s_u = [-sin(lon), -cos(lon)*sin(lat)]
    c2s_v = [ cos(lon), -sin(lon)*sin(lat)]
    c2s_w = [0.       ,           cos(lat)]
    
    us = c2s_u(1)*uc + c2s_v(1)*vc + c2s_w(1)*wc
    vs = c2s_u(2)*uc + c2s_v(2)*vc + c2s_w(2)*wc
  
  end subroutine wind_cart2sph
  
  function point2plane(p1,p2,p3,pd)
    real(r8),dimension(3)             :: point2plane
    real(r8),dimension(3),intent(in ) :: p1,p2,p3
    real(r8),dimension(3),intent(in ) :: pd
    
    real(r8),dimension(3) :: v1,v2,n
    real(r8) A,B,C,x1,y1,z1,xd,yd,zd
    real(r8) lambda
    
    v1 = p2 - p1
    v2 = p3 - p1
    n  = cross(v1,v2)
    A  = n(1)
    B  = n(2)
    C  = n(3)
    
    x1 = p1(1)
    y1 = p1(2)
    z1 = p1(3)
    
    xd = pd(1)
    yd = pd(2)
    zd = pd(3)
    
    lambda = (A*(x1-xd)+B*(y1-yd)+C*(z1-zd))/(A**2+B**2+C**2)
    
    point2plane(1) = xd + A*lambda
    point2plane(2) = yd + B*lambda
    point2plane(3) = zd + C*lambda

  end function point2plane
  
  function point2normalPlane(n,pd)
    real(r8),dimension(3)             :: point2normalPlane
    real(r8),dimension(3),intent(in ) :: n
    real(r8),dimension(3),intent(in ) :: pd
    
    real(r8) A,B,C,x1,y1,z1,xd,yd,zd
    real(r8) lambda
    
    A  = n(1)
    B  = n(2)
    C  = n(3)
    
    x1 = n(1)
    y1 = n(2)
    z1 = n(3)
    
    xd = pd(1)
    yd = pd(2)
    zd = pd(3)
    
    lambda = (A*(x1-xd)+B*(y1-yd)+C*(z1-zd))/(A**2+B**2+C**2)
    
    point2normalPlane(1) = xd + A*lambda
    point2normalPlane(2) = yd + B*lambda
    point2normalPlane(3) = zd + C*lambda

  end function point2normalPlane
  
  ! spherical distance on unit sphere
  function spherical_distance(lat1,lon1,lat2,lon2,r)
    real(r8) :: spherical_distance
    real(r8),intent(in) :: lat1,lon1,lat2,lon2
    real(r8),intent(in) :: r
    
    spherical_distance = r * acos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon1-lon2) )
  end function spherical_distance
  
  ! Calculate spherical polygon area on unit sphere.
  ! Achieved from Dong Li https://github.com/dongli/fishman/blob/master/src/sphere_geometry_mod.F90#L198
  real(r8) function calc_polygon_area(x, y, z) result(res)

    real(r8), intent(in) :: x(:)
    real(r8), intent(in) :: y(:)
    real(r8), intent(in) :: z(:)

    integer n, im1, i, ip1
    real(r8) angle

    n = size(x)
    if (n < 3) then
      print*,'Spherical polygon number is less than 3!'
    end if
    
    res = 0._r16
    do i = 1, n
      im1 = merge(i - 1, n, i /= 1)
      ip1 = merge(i + 1, 1, i /= n)
      angle = calc_sphere_angle([x(im1),y(im1),z(im1)], [x(i),y(i),z(i)], [x(ip1),y(ip1),z(ip1)])
      res = res + angle
    end do
    res = res - ( n - 2 ) * pi
    
    if (abs(res) > 1._r16 .or. res <= 0._r16) then
      print*, 'Encounter bad spherical polygon to calculate area!'
    end if
    
  end function calc_polygon_area
  
  ! Calculate the dihedra angle between plane AB and plane BC.
  ! Achieved from Dong Li https://github.com/dongli/fishman/blob/master/src/sphere_geometry_mod.F90#L247
  real(r8) function calc_sphere_angle(a, b, c) result(res)

    real(r8), intent(in) :: a(3)
    real(r8), intent(in) :: b(3)
    real(r8), intent(in) :: c(3)

    real(r8) nab(3) ! Normal vector of plane AB
    real(r8) nbc(3) ! Normal vector of plane BC

    nab = unify_vector(cross(a, b))
    nbc = unify_vector(cross(b, c))
    res = acos(- max(min(dot_product(nab, nbc), 1._r16), -1._r16))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    !if( dot_product(cross(nab, nbc), a) < 0 ) res = pi2 - res

  end function calc_sphere_angle
  
  ! Calculate unit vector
  function unify_vector(x) result(res)

    real(r8), intent(in) :: x(:)
    real(r8) res(size(x))

    real(r8) n

    n = norm(x)
    if (n /= 0._r16) then
      res = x / n
    else
      res = x
    end if

  end function unify_vector
  
  ! Judge side of point pd on line p1->p2, based on Dongli TTS-I/Sources/Sphere/Sphere.cpp
  ! https://github.com/dongli/TTS-I/blob/master/Sources/Sphere/Sphere.cpp#L444
  function orient(p1,p2,pd)
    integer (i4 ) :: orient
    real    (r8) :: p1(3),p2(3),pd(3)
    
    real(r8) :: det
    
    det = dot_product(pd,cross(p1,p2))

     if(det>0._r16)then
       orient = orientLeft
     elseif(det<0._r16)then
       orient = orientRight
     else
       orient = orientOn
     endif
    
  end function orient
  
  ! Judge whether the point is in given triangle, based on Dongli TTS-I/Sources/Sphere/Sphere.cpp
  ! https://github.com/dongli/TTS-I/blob/master/Sources/Sphere/Sphere.cpp#L609
  function inTriangle(vertex1,vertex2,vertex3,point)
    logical                 :: inTriangle
    real   (r8),intent(in) :: vertex1(3)
    real   (r8),intent(in) :: vertex2(3)
    real   (r8),intent(in) :: vertex3(3)
    real   (r8),intent(in) :: point  (3)
    
    integer(i4) :: vertexOrder
    integer(i4) :: ort(3)
    
    ! if vertexOrder == orientLeft, anti-clock squence
    ! if vertexOrder == orientRight, clock squence
    ! if vertexOrder == orientOn, triangle is too small
    vertexOrder = orient(vertex2,vertex3,vertex1)
    
    if(vertexOrder == orientOn) stop 'triangle is too small'
    
    ort(1) = orient(vertex1,vertex2,point)
    ort(2) = orient(vertex2,vertex3,point)
    ort(3) = orient(vertex3,vertex1,point)
  
    if( (any(ort==orientRight) .and. vertexOrder==orientLeft )&
    .or.(any(ort==orientLeft ) .and. vertexOrder==orientRight) )then
      inTriangle = .false.
    else
      inTriangle = .true.
    endif
    
  end function inTriangle
  
  ! calculate determinant
  function det(matrix)
    real(r8)            :: det
    real(r8),intent(in) :: matrix(:,:)
    
    integer ids,ide
    integer order
    
    ids = lbound(matrix,1)
    ide = ubound(matrix,1)
    
    order = ide - ids + 1
    
    if(order == 2)then
      det = matrix(1,1) * matrix(2,2) - matrix(2,1) * matrix(1,2)
    elseif(order == 3)then
    det = matrix(1,1)*matrix(2,2)*matrix(3,3) + matrix(1,2)*matrix(2,3)*matrix(3,1) + matrix(1,3)*matrix(2,1)*matrix(3,2)&
        - matrix(1,3)*matrix(2,2)*matrix(3,1) - matrix(1,2)*matrix(2,1)*matrix(3,3) - matrix(1,1)*matrix(2,3)*matrix(3,2)
    else
      stop 'Function det do not support determinant over 3nd order'
    endif
  end function det
  
  ! calculate determinant
  recursive function determinant(matrix) result(laplace_det)
    real   (r8)              :: matrix(:,:)
    integer(i4 )              :: i, n, p, q
    real   (r8)              :: laplace_det, det
    real   (r8), allocatable :: cf(:,:)
    
    n = size(matrix, 1) 
    
    !!***** output *****   
    !print "(a)", "Entering determinant() with matrix = "
    !do p = 1, n
    !    print "(4x,100(f3.1,x))", ( matrix( p, q ), q=1,n )
    !enddo
    
    if (n == 1) then  
      det = matrix(1,1)
    else
      det = 0
      do i = 1, n  
        allocate( cf(n-1, n-1) )
        cf = cofactor( matrix, i, 1 )
        
        !!***** output *****
        !print "(4x,a,i0,a,i0,a)", "Getting a ", &
        !        n-1, "-by-", n-1, " sub-matrix from cofactor():"
        !do p = 1, n-1
        !    print "(8x, 100(f3.1,x))", ( cf( p, q ), q=1,n-1 )
        !enddo
        !
        !print "(4x,a)", "and passing it to determinant()."
        
        det = det + ((-1)**(i+1))* matrix( i, 1 ) * determinant( cf )
        deallocate(cf)
      end do
    end if
    
    laplace_det = det
    
    !!***** output *****
    !print *, "  ---> Returning det = ", det
  end function
  
  function cofactor(matrix, mI, mJ)
    real   (r8), dimension(:,:)              :: matrix
    integer(i4 )                              :: mI, mJ
    integer(i4 )                              :: msize(2), i, j, k, l, n
    real   (r8), dimension(:,:), allocatable :: cofactor
    
    msize = shape(matrix)
    n = msize(1)
  
    allocate(cofactor(n-1, n-1))
    l=0
    k = 1
    do i=1, n
     if (i .ne. mI) then
       l = 1
       do j=1, n
         if (j .ne. mJ) then
           cofactor(k,l) = matrix(i,j)
           l = l+ 1
         end if
       end do
       k = k+ 1
     end if
    end do
  end function cofactor
  
  function center_difference(f,dh)
    real              :: center_difference
    real,dimension(:) :: f
    real              :: dh
    
    integer ids,ide
    integer ic
    
    ids = lbound(f,1)
    ide = ubound(f,1)
    
    ic = ( ide - ids + 1 ) / 2 + 1
    
    !! 2nd
    !center_difference = ( f(ic+1) - f(ic-1) ) / ( 2. * dh )
    !!print*,'2nd',center_difference
    
    !! 4th
    !center_difference = ( f(ic-2) - 8.*f(ic-1) + 8.*f(ic+1) - f(ic+2) ) / ( 12. * dh )
    !!print*,'4th',center_difference
    
    ! 6th
    center_difference = ( -f(ic-3) + 9. * f(ic-2) - 45. * f(ic-1) + 45. * f(ic+1) - 9. * f(ic+2) + f(ic+3) ) / ( 60. * dh );
    !print*,'6th',center_difference
    
    !print*,''
    
  end function center_difference
  
      ! calculate inverse matrix of A_input
      ! N is the order of matrix A_input and A
      ! A is inverse A_input
      ! L is status info
	  SUBROUTINE BRINV(N,A_input,A,L)
      implicit none
      integer(i4 ),intent(in ) :: N
	  real   (r8),intent(in ) :: A_input(N,N)
	  real   (r8),intent(out) :: A      (N,N)
      integer(i4 ),intent(out) :: L
      
	  real   (r8) :: T,D
      integer(i4 ) :: IS(N),JS(N)
      integer(i4 ) :: i,j,k
      
      A = A_input
      
	  L=1
	  DO 100 K=1,N
	    D=0._r16
	    DO 10 I=K,N
	    DO 10 J=K,N
	      IF (ABS(A(I,J)).GT.D) THEN
	        D=ABS(A(I,J))
	        IS(K)=I
	        JS(K)=J
	      END IF
10	    CONTINUE
	    IF (D+1.0.EQ.1.0) THEN
	      L=0
	      WRITE(*,20)
	      RETURN
	    END IF
20	    FORMAT(1X,'ERR**NOT INV')
	    DO 30 J=1,N
	      T=A(K,J)
	      A(K,J)=A(IS(K),J)
	      A(IS(K),J)=T
30	    CONTINUE
	    DO 40 I=1,N
	      T=A(I,K)
	      A(I,K)=A(I,JS(K))
	      A(I,JS(K))=T
40	    CONTINUE
	    A(K,K)=1/A(K,K)
	    DO 50 J=1,N
	      IF (J.NE.K) THEN
	        A(K,J)=A(K,J)*A(K,K)
	      END IF
50	    CONTINUE
	    DO 70 I=1,N
	      IF (I.NE.K) THEN
	        DO 60 J=1,N
	          IF (J.NE.K) THEN
	            A(I,J)=A(I,J)-A(I,K)*A(K,J)
	          END IF
60	        CONTINUE
	      END IF
70	    CONTINUE
	    DO 80 I=1,N
	      IF (I.NE.K) THEN
	        A(I,K)=-A(I,K)*A(K,K)
	      END IF
80	    CONTINUE
100	  CONTINUE
	  DO 130 K=N,1,-1
	    DO 110 J=1,N
	      T=A(K,J)
	      A(K,J)=A(JS(K),J)
	      A(JS(K),J)=T
110	    CONTINUE
	    DO 120 I=1,N
	      T=A(I,K)
	      A(I,K)=A(I,IS(K))
	      A(I,IS(K))=T
120	    CONTINUE
130	  CONTINUE
	  RETURN
      END SUBROUTINE BRINV
  
end module math_mod
    