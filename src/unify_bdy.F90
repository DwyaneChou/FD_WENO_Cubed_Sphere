MODULE unify_bdy_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use projection_mod
  
  implicit none
  
    contains
    
    subroutine unify_bdy_field_scalar(field)
      real,intent(inout) :: field(ids:ide,jds:jde,ifs:ife)
      
      real field_raw(ids:ide,jds:jde,ifs:ife)
      
      field_raw = field
      
      ! Bdy lines
      ! edge(1,4)
      field(ids,jds:jde,1) = 0.5 * (field_raw(ide,jds:jde,4) + field_raw(ids,jds:jde,1)) ! left bdy of patch 1
      field(ide,jds:jde,4) = field(ids,jds:jde,1)
      
      ! edge(1,2)
      field(ids,jds:jde,2) = 0.5 * (field_raw(ide,jds:jde,1) + field_raw(ids,jds:jde,2)) ! left bdy of patch 2
      field(ide,jds:jde,1) = field(ids,jds:jde,2)
      
      ! edge(2,3)
      field(ids,jds:jde,3) = 0.5 * (field_raw(ide,jds:jde,2) + field_raw(ids,jds:jde,3)) ! left bdy of patch 3
      field(ide,jds:jde,2) = field(ids,jds:jde,3)
      
      ! edge(3,4)
      field(ids,jds:jde,4) = 0.5 * (field_raw(ide,jds:jde,3) + field_raw(ids,jds:jde,4)) ! left bdy of patch 4
      field(ide,jds:jde,3) = field(ids,jds:jde,4)
      
      ! edge(1,5)
      field(ids:ide,jde,1) = 0.5 * (field_raw(ids:ide,jds,5) + field_raw(ids:ide,jde,1)) ! top bdy of patch 1
      field(ids:ide,jds,5) = field(ids:ide,jde,1)
      
      ! edge(1,6)
      field(ids:ide,jds,1) = 0.5 * (field_raw(ids:ide,jds,1) + field_raw(ids:ide,jde,6)) ! bottom bdy of patch 1
      field(ids:ide,jde,6) = field(ids:ide,jds,1)
      
      ! edge(2,5)
      field(ids:ide,jde,2) = 0.5 * (field_raw(ide,jds:jde,5) + field_raw(ids:ide,jde,2)) ! top bdy of patch 2
      field(ide,jds:jde,5) = field(ids:ide,jde,2)
      
      ! edge(2,6)
      field(ids:ide,jds   ,2) = 0.5 * (field_raw(ids:ide,jds,2) + field_raw(ide,jde:jds:-1,6)) ! bottom bdy of patch 2
      field(ide,jde:jds:-1,6) = field(ids:ide,jds,2)
      
      ! edge(3,5)
      field(ids:ide,jde   ,3) = 0.5 * (field_raw(ids:ide,jde,3) + field_raw(ide:ids:-1,jde,5)) ! top bdy of patch 3
      field(ide:ids:-1,jde,5) = field(ids:ide,jde,3)
      
      ! edge(3,6)
      field(ids:ide,jds   ,3) = 0.5 * (field_raw(ids:ide,jds,3) + field_raw(ide:ids:-1,jds,6)) ! bottom bdy of patch 3
      field(ide:ids:-1,jds,6) = field(ids:ide,jds,3)
      
      ! edge(4,5)
      field(ids:ide,jde   ,4) = 0.5 * (field_raw(ids:ide,jde,4) + field_raw(ids,jde:jds:-1,5)) ! top bdy of patch 4
      field(ids,jde:jds:-1,5) = field(ids:ide,jde,4)
      
      ! edge(4,6)
      field(ids:ide,jds,4) = 0.5 * (field_raw(ids:ide,jds,4) + field_raw(ids,jds:jde,6)) ! bottom bdy of patch 4
      field(ids,jds:jde,6) = field(ids:ide,jds,4)
      
    end subroutine unify_bdy_field_scalar
    
    subroutine unify_bdy_field_contravariant(field_x,field_y)
      real,intent(inout) :: field_x(ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: field_y(ids:ide,jds:jde,ifs:ife)
      
      real us (ids:ide),vs (ids:ide)
      real us1(ids:ide),vs1(ids:ide)
      real us2(ids:ide),vs2(ids:ide)
      
      integer i,j,iPatch
      integer ir,jr
      integer iPatch1, iPatch2
      
      real field_x_raw(ids:ide,jds:jde,ifs:ife)
      real field_y_raw(ids:ide,jds:jde,ifs:ife)
      
      field_x_raw = field_x
      field_y_raw = field_y
      
      ! Bdy lines
      ! edge(1,4)
      iPatch1 = 1
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 4
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo

      ! edge(1,2)
      iPatch1 = 2
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,3)
      iPatch1 = 3
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 2
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,4)
      iPatch1 = 4
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 3
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,5)
      iPatch1 = 5
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,6)
      iPatch1 = 1
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,5)
      iPatch1 = 2
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,6)
      iPatch1 = 2
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        jr = jde - j + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,5)
      iPatch1 = 3
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        ir = ide - i + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,6)
      iPatch1 = 3
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jds
      do i = ids, ide
        ir = ide - i + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,5)
      iPatch1 = 4
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        jr = jde - j + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,6)
      iPatch1 = 4
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
    end subroutine unify_bdy_field_contravariant
    
    subroutine unify_bdy_field_covariant(field_x,field_y)
      real,intent(inout) :: field_x(ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: field_y(ids:ide,jds:jde,ifs:ife)
      
      real us (ids:ide),vs (ids:ide)
      real us1(ids:ide),vs1(ids:ide)
      real us2(ids:ide),vs2(ids:ide)
      
      integer i,j,iPatch
      integer ir,jr
      integer iPatch1, iPatch2
      
      real field_x_raw(ids:ide,jds:jde,ifs:ife)
      real field_y_raw(ids:ide,jds:jde,ifs:ife)
      
      field_x_raw = field_x
      field_y_raw = field_y
      
      ! Bdy lines
      ! edge(1,4)
      iPatch1 = 1
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 4
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      

      ! edge(1,2)
      iPatch1 = 2
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,3)
      iPatch1 = 1
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 4
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,4)
      iPatch1 = 4
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 3
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,5)
      iPatch1 = 5
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,6)
      iPatch1 = 1
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,5)
      iPatch1 = 2
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,6)
      iPatch1 = 2
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        jr = jde - j + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,5)
      iPatch1 = 3
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        ir = ide - i + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,6)
      iPatch1 = 3
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jds
      do i = ids, ide
        ir = ide - i + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,5)
      iPatch1 = 4
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        jr = jde - j + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,6)
      iPatch1 = 4
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
    end subroutine unify_bdy_field_covariant
    
    
END MODULE unify_bdy_mod

