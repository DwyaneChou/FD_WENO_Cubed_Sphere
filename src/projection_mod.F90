MODULE projection_mod
  use constants_mod
  use parameters_mod
  implicit none

contains
  
  subroutine covProjSphere2Plane(cov1, cov2, sv1, sv2, matrixIA, matrixG)

    implicit none
    
    real, intent(out) :: cov1,cov2
    real, intent(in ) :: sv1,sv2
    real, intent(in ) :: matrixIA(2,2)
    real, intent(in ) :: matrixG (2,2)
    
    real matrix(2,2)
    
    matrix = matmul(matrixG,matrixIA)
          
    cov1 = matrix(1,1)*sv1 + matrix(1,2)*sv2
    cov2 = matrix(2,1)*sv1 + matrix(2,2)*sv2

    return
  end subroutine covProjSphere2Plane
  
  subroutine covProjPlane2Sphere(sv1, sv2, cov1, cov2, matrixA, matrixIG)

    implicit none

    real, intent(out) :: sv1,sv2
    real, intent(in)  :: cov1,cov2
    real, intent(in)  :: matrixA (2,2)
    real, intent(in)  :: matrixIG(2,2)
    
    real matrix(2,2)
    
    matrix = matmul(matrixA,matrixIG)
    
    sv1 = matrix(1,1) * cov1 + matrix(1,2) * cov2
    sv2 = matrix(2,1) * cov1 + matrix(2,2) * cov2

    return
  end subroutine covProjPlane2Sphere
  
  subroutine contravProjSphere2Plane(contrav1, contrav2, sv1, sv2, matrixIA)

    implicit none
    
    real   ,intent(in)  :: sv1,sv2
    real   ,intent(in)  :: matrixIA(2,2)
    real   ,intent(out) :: contrav1,contrav2
          
    contrav1 = matrixIA(1,1)*sv1 + matrixIA(1,2)*sv2
    contrav2 = matrixIA(2,1)*sv1 + matrixIA(2,2)*sv2

    return
  end subroutine contravProjSphere2Plane

  subroutine contravProjPlane2Sphere(sv1, sv2, contrav1, contrav2, matrixA)

    implicit none

    real,intent(in)  :: contrav1,contrav2
    real,intent(in)  :: matrixA(2,2)
    real,intent(out) :: sv1,sv2
    
    sv1 = matrixA(1,1) * contrav1 + matrixA(1,2) * contrav2
    sv2 = matrixA(2,1) * contrav1 + matrixA(2,2) * contrav2

    return
  end subroutine contravProjPlane2Sphere

  subroutine contrav2cov(cov1,cov2,contrav1,contrav2, matrixG)
    real, intent(in ) :: contrav1
    real, intent(in ) :: contrav2
    real, intent(in ) :: matrixG(2,2)
    real, intent(out) :: cov1
    real, intent(out) :: cov2
    
    cov1 = matrixG(1,1) * contrav1 + matrixG(1,2) * contrav2
    cov2 = matrixG(2,1) * contrav1 + matrixG(2,2) * contrav2
    
  end subroutine contrav2cov
  
  subroutine cov2contrav(contrav1,contrav2,cov1,cov2,matrixIG)
    real, intent(in ) :: cov1
    real, intent(in ) :: cov2
    real, intent(in ) :: matrixIG(2,2)
    real, intent(out) :: contrav1
    real, intent(out) :: contrav2
    
    contrav1 = matrixIG(1,1) * cov1 + matrixIG(1,2) * cov2
    contrav2 = matrixIG(2,1) * cov1 + matrixIG(2,2) * cov2
    
  end subroutine cov2contrav
  
END MODULE projection_mod

