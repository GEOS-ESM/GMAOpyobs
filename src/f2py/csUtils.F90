module Cubed_Sphere_Utilities
   use, intrinsic :: iso_fortran_env, only: REAL64,REAL32

implicit none
private
public get_cubed_sphere_index
public get_cubed_sphere_index_cart
public get_cubed_sphere_index_convert_to_cart
contains

 subroutine get_cubed_sphere_index_convert_to_cart(center_lons,center_lats,corner_lons,corner_lats,lons,lats,indices,rc)
    real(real64), intent(in) :: center_lats(:,:,:),center_lons(:,:,:)
    real(real64), intent(in) :: corner_lats(:,:,:),corner_lons(:,:,:)
    real(real64), intent(in) :: lons(:),lats(:)
    integer, intent(out) :: indices(:,:)
    integer, intent(out), optional :: rc

    real(real64), allocatable :: centers(:,:,:,:)
    real(real64), allocatable :: corners(:,:,:,:)
    integer :: n,nd,i,j,nf,ii(1),jj(1),im,jm
    logical :: in_region
    real(real64) :: p1(3),p2(3),p3(3),p4(3),p0(3)

    nd=size(corner_lons,3)
    im=size(center_lons,1)
    jm=size(center_lons,2)
    allocate(centers(3,im,jm,nd))
    allocate(corners(3,im+1,jm+1,nd))
    
    do n=1,nd
       do i=1,im
          do j=1,jm
             centers(:,i,j,n)=convert_to_cart([center_lons(i,j,n),center_lats(i,j,n)])
          enddo
       enddo
    enddo
    do n=1,nd
       do i=1,im+1
          do j=1,jm+1
             corners(:,i,j,n)=convert_to_cart([corner_lons(i,j,n),corner_lats(i,j,n)])
          enddo
       enddo
    enddo


    
    do i=1,size(lons)
       do n=1,nd
          p1=corners(:,1,1,n)
          p2=corners(:,im+1,1,n)
          p3=corners(:,im+1,jm+1,n)
          p4=corners(:,1,jm+1,n)
          p0=convert_to_cart([lons(i),lats(i)])
          in_region=point_in_polygon_cart(p0,centers(:,1,1,n),&
               p1,p2,p3,p4)
          if (in_region) then
             nf=n
             exit
          end if
       enddo
       call get_points_in_spherical_domain_cart(centers(:,:,:,nf),corners(:,:,:,nf),[lons(i)],[lats(i)],ii,jj,rc)
       indices(i,1)=ii(1)
       indices(i,2)=jj(1)
       indices(i,3)=nf
    enddo
       
 end subroutine get_cubed_sphere_index_convert_to_cart

 subroutine get_cubed_sphere_index_cart(centers,corners,lons,lats,indices,rc)
    real(real64), intent(in) :: centers(:,:,:,:)
    real(real64), intent(in) :: corners(:,:,:,:)
    real(real64), intent(in) :: lons(:),lats(:)
    integer, intent(out) :: indices(:,:)
    integer, intent(out), optional :: rc

    integer :: n,nd,i,nf,ii(1),jj(1),im,jm
    logical :: in_region
    real(real64) :: p1(3),p2(3),p3(3),p4(3),p0(3)

    nd=size(corners,4)
    im=size(centers,2)
    jm=size(centers,3)
    
    do i=1,size(lons)
       do n=1,nd
          p1=corners(:,1,1,n)
          p2=corners(:,im+1,1,n)
          p3=corners(:,im+1,jm+1,n)
          p4=corners(:,1,jm+1,n)
          p0=convert_to_cart([lons(i),lats(i)])
          in_region=point_in_polygon_cart(p0,centers(:,1,1,n),&
               p1,p2,p3,p4)
          if (in_region) then
             nf=n
             exit
          end if
       enddo
       call get_points_in_spherical_domain_cart(centers(:,:,:,nf),corners(:,:,:,nf),[lons(i)],[lats(i)],ii,jj,rc)
       indices(i,1)=ii(1)
       indices(i,2)=jj(1)
       indices(i,3)=nf
    enddo
       
 end subroutine get_cubed_sphere_index_cart

 subroutine get_cubed_sphere_index(center_lons,center_lats,corner_lons,corner_lats,lons,lats,indices,rc)
    real(real64), intent(in) :: center_lats(:,:,:),center_lons(:,:,:)
    real(real64), intent(in) :: corner_lats(:,:,:),corner_lons(:,:,:)
    real(real64), intent(in) :: lons(:),lats(:)
    integer, intent(out) :: indices(:,:)
    integer, intent(out), optional :: rc

    integer :: n,nd,i,nf,ii(1),jj(1),im,jm
    logical :: in_region
    real(real64) :: p1(2),p2(2),p3(2),p4(2)

    nd=size(corner_lons,3)
    im=size(center_lons,1)
    jm=size(center_lons,2)
    
    do i=1,size(lons)
       do n=1,nd
          p1=[corner_lons(1,1,n),corner_lats(1,1,n)]
          p2=[corner_lons(im+1,1,n),corner_lats(im+1,1,n)]
          p3=[corner_lons(im+1,jm+1,n),corner_lats(im+1,jm+1,n)]
          p4=[corner_lons(1,jm+1,n),corner_lats(1,jm+1,n)]
          in_region=point_in_polygon([lons(i),lats(i)],[center_lons(1,1,n),center_lats(1,1,n)],&
               p1,p2,p3,p4)
          if (in_region) then
             nf=n
             exit
          end if
       enddo
       call get_points_in_spherical_domain(center_lons(:,:,nf),center_lats(:,:,nf), &
             corner_lons(:,:,nf),corner_lats(:,:,nf),[lons(i)],[lats(i)],ii,jj,rc)
       indices(i,1)=ii(1)
       indices(i,2)=jj(1)
       indices(i,3)=nf
    enddo
       
 end subroutine get_cubed_sphere_index

 subroutine get_points_in_spherical_domain_cart(centers,corners,lons,lats,ii,jj,rc)
    real(real64), intent(in) :: centers(:,:,:)
    real(real64), intent(in) :: corners(:,:,:)
    real(real64), intent(in) :: lons(:),lats(:)
    integer, intent(out) :: ii(:),jj(:)
    integer, intent(out), optional :: rc

    integer :: npts,i,n,niter,im,jm,ilb,jlb,iub,jub,ifound,jfound
    integer :: lold,uold,lnew,unew
    logical :: in_region,in_sub_region
    real(real64) :: p(3)

    npts = size(lats)

    im=size(corners,2)-1
    jm=size(corners,3)-1
    niter = max(im,jm)

    do i=1,npts
       ifound=-1
       jfound=-1
       ilb=1
       iub=im
       jlb=1
       jub=jm
       p=convert_to_cart([lons(i),lats(i)])
       in_region = point_in_polygon_cart(p,centers(:,ilb,jlb),  &
          corners(:,ilb,jlb), corners(:,iub+1,jlb),corners(:,iub+1,jub+1),corners(:,ilb,jub+1))
       if (in_region) then
          ! bisect first dimension
          lnew=ilb
          unew=iub 
          do n = 1,niter
             lold=lnew
             uold=unew
             unew=lold+(uold-lold)/2
             in_sub_region = point_in_polygon_cart(p,centers(:,lnew,jlb), &
                 corners(:,lnew,jlb),  &
                 corners(:,unew+1,jlb), &
                 corners(:,unew+1,jub+1), &
                 corners(:,lnew,jub+1))
             if (in_sub_region) then
               lnew=lold
               unew=unew 
             else
               lnew=unew+1
               unew=uold
             end if
             if (unew==lnew) then
                ifound=unew
                exit
             end if
          enddo
          ! bisect 2nd dimension
          lnew=jlb
          unew=jub
          do n = 1,niter
             lold=lnew
             uold=unew
             unew=lold+(uold-lold)/2
             in_sub_region = point_in_polygon_cart(p, centers(:,ifound,lnew) , &
                 corners(:,ifound,lnew), &
                 corners(:,ifound+1,lnew),  &
                 corners(:,ifound+1,unew+1), &
                 corners(:,ifound,unew+1))
             if (in_sub_region) then
               lnew=lold
               unew=unew 
             else
               lnew=unew+1
               unew=uold
             end if
             if (unew==lnew) then
                jfound=unew
                exit
             end if
          enddo
       end if
       ii(i)=ifound
       jj(i)=jfound
    enddo
         
 end subroutine get_points_in_spherical_domain_cart 

 subroutine get_points_in_spherical_domain(center_lons,center_lats,corner_lons,corner_lats,lons,lats,ii,jj,rc)
    real(real64), intent(in) :: center_lats(:,:),center_lons(:,:)
    real(real64), intent(in) :: corner_lats(:,:),corner_lons(:,:)
    real(real64), intent(in) :: lons(:),lats(:)
    integer, intent(out) :: ii(:),jj(:)
    integer, intent(out), optional :: rc

    integer :: npts,i,n,niter,im,jm,ilb,jlb,iub,jub,ifound,jfound
    integer :: lold,uold,lnew,unew
    logical :: in_region,in_sub_region

    npts = size(lats)

    im=size(corner_lons,1)-1
    jm=size(corner_lons,2)-1
    niter = max(im,jm)

    do i=1,npts
       ifound=-1
       jfound=-1
       ilb=1
       iub=im
       jlb=1
       jub=jm
       in_region = point_in_polygon([lons(i),lats(i)],[center_lons(ilb,jlb),center_lats(ilb,jlb)],  &
           [corner_lons(ilb,jlb),corner_lats(ilb,jlb)], &
           [corner_lons(iub+1,jlb),corner_lats(iub+1,jlb)], &
           [corner_lons(iub+1,jub+1),corner_lats(iub+1,jub+1)], &
           [corner_lons(ilb,jub+1),corner_lats(ilb,jub+1)])
       if (in_region) then
          ! bisect first dimension
          lnew=ilb
          unew=iub 
          do n = 1,niter
             lold=lnew
             uold=unew
             unew=lold+(uold-lold)/2
             in_sub_region = point_in_polygon([lons(i),lats(i)], [center_lons(lnew,jlb),center_lats(lnew,jlb)], &
                 [corner_lons(lnew,jlb),corner_lats(lnew,jlb)], &
                 [corner_lons(unew+1,jlb),corner_lats(unew+1,jlb)], &
                 [corner_lons(unew+1,jub+1),corner_lats(unew+1,jub+1)], &
                 [corner_lons(lnew,jub+1),corner_lats(lnew,jub+1)])
             if (in_sub_region) then
               lnew=lold
               unew=unew 
             else
               lnew=unew+1
               unew=uold
             end if
             if (unew==lnew) then
                ifound=unew
                exit
             end if
          enddo
          ! bisect 2nd dimension
          lnew=jlb
          unew=jub
          do n = 1,niter
             lold=lnew
             uold=unew
             unew=lold+(uold-lold)/2
             in_sub_region = point_in_polygon([lons(i),lats(i)], [center_lons(ifound,lnew),center_lats(ifound,lnew)] , &
                 [corner_lons(ifound,lnew),corner_lats(ifound,lnew)], &
                 [corner_lons(ifound+1,lnew),corner_lats(ifound+1,lnew)], &
                 [corner_lons(ifound+1,unew+1),corner_lats(ifound+1,unew+1)], &
                 [corner_lons(ifound,unew+1),corner_lats(ifound,unew+1)])
             if (in_sub_region) then
               lnew=lold
               unew=unew 
             else
               lnew=unew+1
               unew=uold
             end if
             if (unew==lnew) then
                jfound=unew
                exit
             end if
          enddo
       end if
       ii(i)=ifound
       jj(i)=jfound
    enddo
         
 end subroutine get_points_in_spherical_domain 

 function point_in_polygon(p0,pinside,a1,a2,a3,a4) result(in_poly)
    real(real64), intent(in) :: p0(2),pinside(2),a1(2),a2(2),a3(2),a4(2)
    logical :: in_poly
 
    real(real64) :: p1c(3),p2c(3),a1c(3),a2c(3),a3c(3),a4c(3)
    logical :: intersect(4)
    p1c=convert_to_cart(p0)
    p2c=convert_to_cart(pinside)
    a1c=convert_to_cart(a1)
    a2c=convert_to_cart(a2)
    a3c=convert_to_cart(a3)
    a4c=convert_to_cart(a4)

    intersect(1) = lines_intersect(p1c,p2c,a1c,a2c)
    intersect(2) = lines_intersect(p1c,p2c,a2c,a3c)
    intersect(3) = lines_intersect(p1c,p2c,a3c,a4c)
    intersect(4) = lines_intersect(p1c,p2c,a4c,a1c)
    if (mod(count(intersect),2)==0) then
       in_poly=.true.
    else
       in_poly=.false.
    end if

 end function point_in_polygon

 function point_in_polygon_cart(p0,p1,a1,a2,a3,a4) result(in_poly)
    real(real64), intent(in) :: p0(3),p1(3),a1(3),a2(3),a3(3),a4(3)
    logical :: in_poly
 
    logical :: intersect(4)

    intersect(1) = lines_intersect(p0,p1,a1,a2)
    intersect(2) = lines_intersect(p0,p1,a2,a3)
    intersect(3) = lines_intersect(p0,p1,a3,a4)
    intersect(4) = lines_intersect(p0,p1,a4,a1)
    if (mod(count(intersect),2)==0) then
       in_poly=.true.
    else
       in_poly=.false.
    end if

 end function point_in_polygon_cart

! it is claimed this should work but doesn't
 !function point_in_polygon_crosprod(p1,a1,a2,a3,a4) result(in_poly)
    !real(real64), intent(in) :: p1(2),a1(2),a2(2),a3(2),a4(2)
    !logical :: in_poly
 
    !real(real64) :: p1c(3),a1c(3),a2c(3),a3c(3),a4c(3)
    !real(real64) :: crs12(3),crs23(3),crs34(3),crs41(3)
    !real(real64) :: d12,d23,d34,d41
    !logical :: signs(4)
    !! a1 -> a2 -> a3 -> a4 so a4 connects to a1

    !p1c=convert_to_cart(p1)
    !a1c=convert_to_cart(a1)
    !a2c=convert_to_cart(a2)
    !a3c=convert_to_cart(a3)
    !a4c=convert_to_cart(a4)

    !crs12 = cross_prod(a1c,a2c)
    !crs23 = cross_prod(a2c,a3c)
    !crs34 = cross_prod(a3c,a4c)
    !crs41 = cross_prod(a4c,a1c)
    !d12=dot_product(p1c,crs12)
    !d23=dot_product(p1c,crs23)
    !d34=dot_product(p1c,crs34)
    !d41=dot_product(p1c,crs41)
    !signs(1)= (d12<0.0)
    !signs(2)= (d23<0.0)
    !signs(3)= (d34<0.0)
    !signs(4)= (d41<0.0)
    !in_poly=( (count(signs)==0) .or. (count(signs)==4) )

 !end function point_in_polygon_crossprod

 function lines_intersect(b0,b1,a0,a1)  result(intersect)
    real(real64), intent(in) :: b0(3),b1(3),a0(3),a1(3)
    logical :: intersect
    real(real64) :: p(3),q(3),t(3)
    real(real64) :: s1,s2,s3,s4
    logical :: signs(4)

    intersect=.false.
    q=cross_prod(b0,b1)
    p=cross_prod(a0,a1)
    t=normal_vect(cross_prod(p,q))

    s1=dot_product(cross_prod(a0,p),t)
    s2=dot_product(cross_prod(a1,p),t)
    s3=dot_product(cross_prod(b0,q),t)
    s4=dot_product(cross_prod(b1,q),t)

    signs(1) = -s1 <0.d0
    signs(2) = s2 <0.d0
    signs(3) = -s3 < 0.d0
    signs(4) = s4 < 0.d0

    intersect = ((count(signs)==0) .or. (count(signs)==4))

 end function lines_intersect

 function normal_vect(vin) result(vout)
    real(real64), intent(in) :: vin(3)
    real(real64) :: vout(3)
    vout=vin/sqrt(vin(1)*vin(1)+vin(2)*vin(2)+vin(3)*vin(3))

 end function normal_vect

 function cross_prod(v1,v2) result(vout)
    real(real64), intent(in) :: v1(3),v2(3)
    real(real64) :: vout(3)
    vout(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vout(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vout(3)=v1(1)*v2(2)-v1(2)*v2(1)
 end function cross_prod

 function convert_to_cart(v) result(xyz)
    real(real64), intent(in) :: v(2)
    real(real64) :: xyz(3)

    xyz(1)=cos(v(2))*cos(v(1))
    xyz(2)=cos(v(2))*sin(v(1))
    xyz(3)=sin(v(2))

 end function convert_to_cart

function vect_mag(v) result(mag)
   real(real64), intent(in) :: v(3)
   real :: mag
   mag = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
end function vect_mag

end module Cubed_Sphere_Utilities
