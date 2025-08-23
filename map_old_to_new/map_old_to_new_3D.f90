program map_old_to_new_3d_nonuniform
    implicit none

    !--------------------------- Variable Declarations ----------------------------
    character(len=100), parameter :: oldFile = "20u.bin" ! Assumed to have u,v,w,T
    character(len=100), parameter :: newFile = "1134.bin"

    ! Grid dimensions
    integer, parameter :: nx_old = 21, ny_old = nx_old, nz_old = nx_old
    integer, parameter :: nx_new = 41, ny_new = nx_new, nz_new = nx_new

    ! Grid distribution types
    logical, parameter :: old_uniform = .true.  ! Set to .true. for uniform grid
    logical, parameter :: new_uniform = .true.  ! Set to .true. for uniform grid

    ! Constants for grid distribution
    real(8), parameter :: stretch_factor_old = 2.1d0
    real(8), parameter :: stretch_factor_new = 2.1d0

    ! Coordinate arrays for grids
    real(8), allocatable :: x_old(:), y_old(:), z_old(:)
    real(8), allocatable :: x_new(:), y_new(:), z_new(:)

    ! Data arrays for old and new grids
    real(8), allocatable :: u_old(:,:,:), v_old(:,:,:), w_old(:,:,:), T_old(:,:,:)
    real(8), allocatable :: u_new(:,:,:), v_new(:,:,:), w_new(:,:,:), T_new(:,:,:)
    integer :: ios, i, j, k

    integer, parameter :: read_unit = 10, write_unit = 20

    character(len=100) :: baseName,oldPltFile, newPltFile

    !--------------------------- Generate Mesh ----------------------------
    allocate(x_old(nx_old), y_old(ny_old), z_old(nz_old))
    allocate(x_new(nx_new), y_new(ny_new), z_new(nz_new))

    ! Compute nonuniform grid coordinates for old grid
    baseName = TRIM(oldFile)
    oldPltFile = baseName(1:INDEX(baseName, '.bin')-1) // '_old.plt'
    baseName = TRIM(newFile)
    newPltFile = baseName(1:INDEX(baseName, '.bin')-1) // '_new.plt'    
    
    write(*,*) "old Tecplot output: ", trim(oldPltFile)
    write(*,*) "new Tecplot output: ", trim(newPltFile)

    if (old_uniform .eqv. .true.) then
        do i = 1, nx_old
            x_old(i) = (dble(i) - 0.5d0) / dble(nx_old)
        end do
        do j = 1, ny_old
            y_old(j) = (dble(j) - 0.5d0) / dble(ny_old)
        end do
        do k = 1, nz_old
            z_old(k) = (dble(k) - 0.5d0) / dble(nz_old)
        end do
    else
        do i = 1, nx_old
            x_old(i) = 0.5d0 * (erf(stretch_factor_old * (dble(i) / dble(nx_old+1) - 0.5d0)) / &
                             erf(0.5d0 * stretch_factor_old) + 1.0d0)
        end do
        do j = 1, ny_old
            y_old(j) = 0.5d0 * (erf(stretch_factor_old * (dble(j) / dble(ny_old+1) - 0.5d0)) / &
                             erf(0.5d0 * stretch_factor_old) + 1.0d0)
        end do
        do k = 1, nz_old
            z_old(k) = 0.5d0 * (erf(stretch_factor_old * (dble(k) / dble(nz_old+1) - 0.5d0)) / &
                             erf(0.5d0 * stretch_factor_old) + 1.0d0)
        end do
    end if

    ! Compute grid coordinates for new grid
    if (new_uniform .eqv. .true.) then
        do i = 1, nx_new
            x_new(i) = (dble(i) - 0.5d0) / dble(nx_new)
        end do
        do j = 1, ny_new
            y_new(j) = (dble(j) - 0.5d0) / dble(ny_new)
        end do
        do k = 1, nz_new
            z_new(k) = (dble(k) - 0.5d0) / dble(nz_new)
        end do
    else
        do i = 1, nx_new
            x_new(i) = 0.5d0 * (erf(stretch_factor_new * (dble(i) / dble(nx_new+1) - 0.5d0)) / &
                             erf(0.5d0 * stretch_factor_new) + 1.0d0)
        end do
        do j = 1, ny_new
            y_new(j) = 0.5d0 * (erf(stretch_factor_new * (dble(j) / dble(ny_new+1) - 0.5d0)) / &
                             erf(0.5d0 * stretch_factor_new) + 1.0d0)
        end do
        do k = 1, nz_new
            z_new(k) = 0.5d0 * (erf(stretch_factor_new * (dble(k) / dble(nz_new+1) - 0.5d0)) / &
                             erf(0.5d0 * stretch_factor_new) + 1.0d0)
        end do
    end if
    !--------------------------- Generate Mesh ----------------------------

    !--------------------------- Data Input ----------------------------
    allocate(u_old(nx_old,ny_old,nz_old), v_old(nx_old,ny_old,nz_old), w_old(nx_old,ny_old,nz_old), &
             T_old(nx_old,ny_old,nz_old))
    open(unit=read_unit, file=oldFile, form="unformatted", access="sequential", status='old', iostat=ios)
    read(read_unit) (((u_old(i,j,k), i=1,nx_old), j=1,ny_old), k=1,nz_old)
    read(read_unit) (((v_old(i,j,k), i=1,nx_old), j=1,ny_old), k=1,nz_old)
    read(read_unit) (((w_old(i,j,k), i=1,nx_old), j=1,ny_old), k=1,nz_old)
    read(read_unit) (((T_old(i,j,k), i=1,nx_old), j=1,ny_old), k=1,nz_old)
    close(read_unit)

    !--------------------------- Check for Errors ----------------------------
    write(*,*) "old data loaded successfully!"
    if (old_uniform .eqv. .false.) then
        write(*,*) "old mesh is nonuniform"
        write(*,*) "stretch factor (old)=", real(stretch_factor_old)
    else
        write(*,*) "old mesh is uniform; NO stretch factor"
    end if
    write(*,*) "nx_old=", nx_old, " ny_old=", ny_old, " nz_old=", nz_old
    write(*,*) "u_old (center)=", real(u_old((nx_old-1)/2+1,(ny_old-1)/2+1,(nz_old-1)/2+1))
    write(*,*) "v_old (center)=", real(v_old((nx_old-1)/2+1,(ny_old-1)/2+1,(nz_old-1)/2+1))
    write(*,*) "w_old (center)=", real(w_old((nx_old-1)/2+1,(ny_old-1)/2+1,(nz_old-1)/2+1))
    write(*,*) "T_old (center)=", real(T_old((nx_old-1)/2+1,(ny_old-1)/2+1,(nz_old-1)/2+1))
    call output_Tecplot_3d(nx_old, ny_old, nz_old, x_old, y_old, z_old, u_old, v_old, w_old, T_old, oldPltFile)
    !--------------------------- Check for Errors ----------------------------

    !--------------------------- Interpolation ----------------------------
    allocate(u_new(nx_new,ny_new,nz_new), v_new(nx_new,ny_new,nz_new), w_new(nx_new,ny_new,nz_new), &
             T_new(nx_new,ny_new,nz_new))
    call conservative_remap_3d( &
        x_old, y_old, z_old, u_old, v_old, w_old, T_old, &
        x_new, y_new, z_new, u_new, v_new, w_new, T_new,&
        nx_old, ny_old, nz_old, nx_new, ny_new, nz_new)

    !--------------------------- Check for Errors ----------------------------
    write(*,*) "Interpolation completed successfully!"
    if (new_uniform .eqv. .false.) then
        write(*,*) "new mesh is nonuniform"
        write(*,*) "stretch factor (new)=", real(stretch_factor_new)
    else
        write(*,*) "new mesh is uniform; NO stretch factor"
    end if
    write(*,*) "nx_new=", nx_new, " ny_new=", ny_new, " nz_new=", nz_new
    write(*,*) "u_new (center)=", real(u_new((nx_new-1)/2+1,(ny_new-1)/2+1,(nz_new-1)/2+1))
    write(*,*) "v_new (center)=", real(v_new((nx_new-1)/2+1,(ny_new-1)/2+1,(nz_new-1)/2+1))
    write(*,*) "w_new (center)=", real(w_new((nx_new-1)/2+1,(ny_new-1)/2+1,(nz_new-1)/2+1))
    write(*,*) "T_new (center)=", real(T_new((nx_new-1)/2+1,(ny_new-1)/2+1,(nz_new-1)/2+1))
    call output_Tecplot_3d(nx_new, ny_new, nz_new, x_new, y_new, z_new, u_new, v_new, w_new, T_new,newPltFile)
    !--------------------------- Check for Errors ----------------------------

    !--------------------------- Data Output ----------------------------
    open(unit=write_unit, file=newFile, form="unformatted", access="sequential", status='replace', iostat=ios)
    write(write_unit) (((u_new(i,j,k), i=1,nx_new), j=1,ny_new), k=1,nz_new)
    write(write_unit) (((v_new(i,j,k), i=1,nx_new), j=1,ny_new), k=1,nz_new)
    write(write_unit) (((w_new(i,j,k), i=1,nx_new), j=1,ny_new), k=1,nz_new)
    write(write_unit) (((T_new(i,j,k), i=1,nx_new), j=1,ny_new), k=1,nz_new)
    close(write_unit)

    write(*,*) "new data written successfully!"

    deallocate(x_old, y_old, z_old, x_new, y_new, z_new)
    deallocate(u_old, v_old, w_old, T_old)
    deallocate(u_new, v_new, w_new, T_new)

    print *, "3D mapping completed successfully!"
end program map_old_to_new_3d_nonuniform

subroutine conservative_remap_3d( &
    x_src, y_src, z_src, u_src, v_src, w_src, T_src, &
    x_tgt, y_tgt, z_tgt, u_tgt, v_tgt, w_tgt, T_tgt, &
    nx_src, ny_src, nz_src, nx_tgt, ny_tgt, nz_tgt)
    
    integer, intent(in) :: nx_src, ny_src, nz_src, nx_tgt, ny_tgt, nz_tgt
    real(8), intent(in) :: x_src(nx_src), y_src(ny_src), z_src(nz_src)
    real(8), intent(in) :: u_src(nx_src, ny_src, nz_src)
    real(8), intent(in) :: v_src(nx_src, ny_src, nz_src)
    real(8), intent(in) :: w_src(nx_src, ny_src, nz_src)
    real(8), intent(in) :: T_src(nx_src, ny_src, nz_src)
    real(8), intent(in) :: x_tgt(nx_tgt), y_tgt(ny_tgt), z_tgt(nz_tgt)
    real(8), intent(out) :: u_tgt(nx_tgt, ny_tgt, nz_tgt)
    real(8), intent(out) :: v_tgt(nx_tgt, ny_tgt, nz_tgt)
    real(8), intent(out) :: w_tgt(nx_tgt, ny_tgt, nz_tgt)
    real(8), intent(out) :: T_tgt(nx_tgt, ny_tgt, nz_tgt)

    ! --- Local variables ---
    integer :: i, j, k, ii, jj, kk
    ! Dynamic search range indices
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    
    real(8) :: total_volume
    real(8) :: u_integral, v_integral, w_integral, T_integral
    real(8) :: overlap_volume, x_overlap_center, y_overlap_center, z_overlap_center
    real(8) :: u_recon, v_recon, w_recon, T_recon

    ! --- Pre-calculated grid properties ---
    real(8), allocatable :: x_src_bounds(:,:), y_src_bounds(:,:), z_src_bounds(:,:)
    real(8), allocatable :: x_tgt_bounds(:,:), y_tgt_bounds(:,:), z_tgt_bounds(:,:)
    
    ! --- Gradients on source grid ---
    real(8), allocatable :: u_grad_x(:,:,:), u_grad_y(:,:,:), u_grad_z(:,:,:)
    real(8), allocatable :: v_grad_x(:,:,:), v_grad_y(:,:,:), v_grad_z(:,:,:)
    real(8), allocatable :: w_grad_x(:,:,:), w_grad_y(:,:,:), w_grad_z(:,:,:)
    real(8), allocatable :: T_grad_x(:,:,:), T_grad_y(:,:,:), T_grad_z(:,:,:)

    ! 1. Allocate space for gradients
    allocate(u_grad_x(nx_src,ny_src,nz_src), u_grad_y(nx_src,ny_src,nz_src), u_grad_z(nx_src,ny_src,nz_src))
    allocate(v_grad_x(nx_src,ny_src,nz_src), v_grad_y(nx_src,ny_src,nz_src), v_grad_z(nx_src,ny_src,nz_src))
    allocate(w_grad_x(nx_src,ny_src,nz_src), w_grad_y(nx_src,ny_src,nz_src), w_grad_z(nx_src,ny_src,nz_src))
    allocate(T_grad_x(nx_src,ny_src,nz_src), T_grad_y(nx_src,ny_src,nz_src), T_grad_z(nx_src,ny_src,nz_src))

    ! 2. Compute limited gradients for all source fields
    call compute_limited_gradients_3d(u_src, x_src, y_src, z_src, nx_src, ny_src, nz_src, u_grad_x, u_grad_y, u_grad_z)
    call compute_limited_gradients_3d(v_src, x_src, y_src, z_src, nx_src, ny_src, nz_src, v_grad_x, v_grad_y, v_grad_z)
    call compute_limited_gradients_3d(w_src, x_src, y_src, z_src, nx_src, ny_src, nz_src, w_grad_x, w_grad_y, w_grad_z)
    call compute_limited_gradients_3d(T_src, x_src, y_src, z_src, nx_src, ny_src, nz_src, T_grad_x, T_grad_y, T_grad_z)

    ! 3. Pre-calculate cell boundaries to improve performance
    allocate(x_src_bounds(2, nx_src), y_src_bounds(2, ny_src), z_src_bounds(2, nz_src))
    allocate(x_tgt_bounds(2, nx_tgt), y_tgt_bounds(2, ny_tgt), z_tgt_bounds(2, nz_tgt))
    call calculate_cell_bounds(x_src, nx_src, x_src_bounds)
    call calculate_cell_bounds(y_src, ny_src, y_src_bounds)
    call calculate_cell_bounds(z_src, nz_src, z_src_bounds)
    call calculate_cell_bounds(x_tgt, nx_tgt, x_tgt_bounds)
    call calculate_cell_bounds(y_tgt, ny_tgt, y_tgt_bounds)
    call calculate_cell_bounds(z_tgt, nz_tgt, z_tgt_bounds)

    ! 4. Main loop over target grid
    do k = 1, nz_tgt
        do j = 1, ny_tgt
            do i = 1, nx_tgt
                i_start = max(1, bounded_bisect(x_src, nx_src, x_tgt_bounds(1, i)) - 1)
                i_end   = min(nx_src, bounded_bisect(x_src, nx_src, x_tgt_bounds(2, i)) + 1)
                j_start = max(1, bounded_bisect(y_src, ny_src, y_tgt_bounds(1, j)) - 1)
                j_end   = min(ny_src, bounded_bisect(y_src, ny_src, y_tgt_bounds(2, j)) + 1)
                k_start = max(1, bounded_bisect(z_src, nz_src, z_tgt_bounds(1, k)) - 1)
                k_end   = min(nz_src, bounded_bisect(z_src, nz_src, z_tgt_bounds(2, k)) + 1)
                
                u_integral = 0.0d0; v_integral = 0.0d0; w_integral = 0.0d0
                T_integral = 0.0d0; 
                total_volume = 0.0d0

                do kk = k_start, k_end
                    do jj = j_start, j_end
                        do ii = i_start, i_end
                            call calculate_overlap_properties_3d( &
                                x_src_bounds(:, ii), y_src_bounds(:, jj), z_src_bounds(:, kk), &
                                x_tgt_bounds(:, i),  y_tgt_bounds(:, j),  z_tgt_bounds(:, k),  &
                                overlap_volume, x_overlap_center, y_overlap_center, z_overlap_center)

                            if (overlap_volume > 1.0d-14) then
                                u_recon = u_src(ii,jj,kk) + u_grad_x(ii,jj,kk) * (x_overlap_center - x_src(ii)) &
                                                       + u_grad_y(ii,jj,kk) * (y_overlap_center - y_src(jj)) &
                                                       + u_grad_z(ii,jj,kk) * (z_overlap_center - z_src(kk))
                                v_recon = v_src(ii,jj,kk) + v_grad_x(ii,jj,kk) * (x_overlap_center - x_src(ii)) &
                                                       + v_grad_y(ii,jj,kk) * (y_overlap_center - y_src(jj)) &
                                                       + v_grad_z(ii,jj,kk) * (z_overlap_center - z_src(kk))
                                w_recon = w_src(ii,jj,kk) + w_grad_x(ii,jj,kk) * (x_overlap_center - x_src(ii)) &
                                                       + w_grad_y(ii,jj,kk) * (y_overlap_center - y_src(jj)) &
                                                       + w_grad_z(ii,jj,kk) * (z_overlap_center - z_src(kk))
                                T_recon = T_src(ii,jj,kk) + T_grad_x(ii,jj,kk) * (x_overlap_center - x_src(ii)) &
                                                       + T_grad_y(ii,jj,kk) * (y_overlap_center - y_src(jj)) &
                                                       + T_grad_z(ii,jj,kk) * (z_overlap_center - z_src(kk))

                                u_integral = u_integral + u_recon * overlap_volume
                                v_integral = v_integral + v_recon * overlap_volume
                                w_integral = w_integral + w_recon * overlap_volume
                                T_integral = T_integral + T_recon * overlap_volume
                                total_volume = total_volume + overlap_volume
                            end if
                        end do
                    end do
                end do

                if (total_volume > 1.0d-12) then
                    u_tgt(i,j,k) = u_integral / total_volume
                    v_tgt(i,j,k) = v_integral / total_volume
                    w_tgt(i,j,k) = w_integral / total_volume
                    T_tgt(i,j,k) = T_integral / total_volume
                else
                    ii = bounded_bisect(x_src, nx_src, x_tgt(i))
                    jj = bounded_bisect(y_src, ny_src, y_tgt(j))
                    kk = bounded_bisect(z_src, nz_src, z_tgt(k))
                    u_tgt(i,j,k) = u_src(ii,jj,kk)
                    v_tgt(i,j,k) = v_src(ii,jj,kk)
                    w_tgt(i,j,k) = w_src(ii,jj,kk)
                    T_tgt(i,j,k) = T_src(ii,jj,kk)
                end if
            end do
        end do
    end do

    ! 5. Deallocate all temporary arrays
    deallocate(u_grad_x, u_grad_y, u_grad_z, v_grad_x, v_grad_y, v_grad_z, w_grad_x, w_grad_y, w_grad_z)
    deallocate(T_grad_x, T_grad_y, T_grad_z)
    deallocate(x_src_bounds, y_src_bounds, z_src_bounds, x_tgt_bounds, y_tgt_bounds, z_tgt_bounds)

contains

    subroutine compute_limited_gradients_3d(field, x, y, z, nx, ny, nz, grad_x, grad_y, grad_z)
        real(8), intent(in) :: field(nx,ny,nz), x(nx), y(ny), z(nz)
        integer, intent(in) :: nx, ny, nz
        real(8), intent(out) :: grad_x(nx,ny,nz), grad_y(nx,ny,nz), grad_z(nx,ny,nz)
        
        integer :: i, j, k
        real(8) :: slope_L, slope_R, slope_B, slope_T, slope_D, slope_U, slope_cen

        ! Compute gradients in x-direction
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (i > 1 .and. i < nx) then
                        slope_L = (field(i,j,k) - field(i-1,j,k)) / (x(i) - x(i-1))
                        slope_R = (field(i+1,j,k) - field(i,j,k)) / (x(i+1) - x(i))
                        slope_cen = (field(i+1,j,k) - field(i-1,j,k)) / (x(i+1) - x(i-1))
                        grad_x(i,j,k) = van_leer_limiter(slope_cen, slope_L, slope_R)
                    else if (i == 1) then
                        grad_x(i,j,k) = (field(2,j,k) - field(1,j,k)) / (x(2) - x(1))
                    else
                        grad_x(i,j,k) = (field(nx,j,k) - field(nx-1,j,k)) / (x(nx) - x(nx-1))
                    end if
                end do
            end do
        end do
        
        ! Compute gradients in y-direction
        do k = 1, nz
            do i = 1, nx
                do j = 1, ny
                    if (j > 1 .and. j < ny) then
                        slope_B = (field(i,j,k) - field(i,j-1,k)) / (y(j) - y(j-1))
                        slope_T = (field(i,j+1,k) - field(i,j,k)) / (y(j+1) - y(j))
                        slope_cen = (field(i,j+1,k) - field(i,j-1,k)) / (y(j+1) - y(j-1))
                        grad_y(i,j,k) = van_leer_limiter(slope_cen, slope_B, slope_T)
                    else if (j == 1) then
                        grad_y(i,j,k) = (field(i,2,k) - field(i,1,k)) / (y(2) - y(1))
                    else
                        grad_y(i,j,k) = (field(i,ny,k) - field(i,ny-1,k)) / (y(ny) - y(ny-1))
                    end if
                end do
            end do
        end do
        
        ! Compute gradients in z-direction
        do j = 1, ny
            do i = 1, nx
                do k = 1, nz
                    if (k > 1 .and. k < nz) then
                        slope_D = (field(i,j,k) - field(i,j,k-1)) / (z(k) - z(k-1))
                        slope_U = (field(i,j,k+1) - field(i,j,k)) / (z(k+1) - z(k))
                        slope_cen = (field(i,j,k+1) - field(i,j,k-1)) / (z(k+1) - z(k-1))
                        grad_z(i,j,k) = van_leer_limiter(slope_cen, slope_D, slope_U)
                    else if (k == 1) then
                        grad_z(i,j,k) = (field(i,j,2) - field(i,j,1)) / (z(2) - z(1))
                    else
                        grad_z(i,j,k) = (field(i,j,nz) - field(i,j,nz-1)) / (z(nz) - z(nz-1))
                    end if
                end do
            end do
        end do
    end subroutine compute_limited_gradients_3d

    function van_leer_limiter(slope_cen, slope_1, slope_2) result(limited_slope)
        real(8), intent(in) :: slope_cen, slope_1, slope_2
        real(8) :: limited_slope, r
        
        if (slope_1 * slope_2 > 0.0d0) then
            r = slope_1 / slope_2
            limited_slope = slope_cen * ( (r + abs(r)) / (1.0d0 + r*r) )
        else
            limited_slope = 0.0d0
        end if
    end function van_leer_limiter

    subroutine calculate_cell_bounds(arr, n, bounds)
        integer, intent(in) :: n
        real(8), intent(in) :: arr(n)
        real(8), intent(out) :: bounds(2, n)
        integer :: i
        do i = 1, n
            if (i == 1) then
                bounds(1, i) = arr(1) - 0.5d0*(arr(2)-arr(1))
            else
                bounds(1, i) = 0.5d0*(arr(i-1) + arr(i))
            end if
            if (i == n) then
                bounds(2, i) = arr(n) + 0.5d0*(arr(n)-arr(n-1))
            else
                bounds(2, i) = 0.5d0*(arr(i) + arr(i+1))
            end if
        end do
    end subroutine calculate_cell_bounds

    subroutine calculate_overlap_properties_3d(b_src_x, b_src_y, b_src_z, b_tgt_x, b_tgt_y, b_tgt_z, volume, xc, yc, zc)
        real(8), intent(in) :: b_src_x(2), b_src_y(2), b_src_z(2)
        real(8), intent(in) :: b_tgt_x(2), b_tgt_y(2), b_tgt_z(2)
        real(8), intent(out) :: volume, xc, yc, zc
        real(8) :: overlap_dx, overlap_dy, overlap_dz, x_min, x_max, y_min, y_max, z_min, z_max

        x_min = max(b_src_x(1), b_tgt_x(1)); x_max = min(b_src_x(2), b_tgt_x(2))
        y_min = max(b_src_y(1), b_tgt_y(1)); y_max = min(b_src_y(2), b_tgt_y(2))
        z_min = max(b_src_z(1), b_tgt_z(1)); z_max = min(b_src_z(2), b_tgt_z(2))
        
        overlap_dx = x_max - x_min
        overlap_dy = y_max - y_min
        overlap_dz = z_max - z_min

        if (overlap_dx > 0.0d0 .and. overlap_dy > 0.0d0 .and. overlap_dz > 0.0d0) then
            volume = overlap_dx * overlap_dy * overlap_dz
            xc = 0.5d0 * (x_min + x_max)
            yc = 0.5d0 * (y_min + y_max)
            zc = 0.5d0 * (z_min + z_max)
        else
            volume = 0.0d0
            xc = 0.0d0; yc = 0.0d0; zc = 0.0d0
        end if
    end subroutine calculate_overlap_properties_3d

    function bounded_bisect(arr, n, x) result(idx)
        integer, intent(in) :: n
        real(8), intent(in) :: arr(n), x
        integer :: idx, lo, hi, mid
        
        if (x <= arr(1)) then; idx = 1; return; end if
        if (x >= arr(n)) then; idx = n; return; end if

        lo = 1; hi = n
        do while (hi - lo > 1); mid = (lo + hi)/2
            if (x >= arr(mid)) then; lo = mid; else; hi = mid; end if
        end do
        idx = lo
    end function bounded_bisect

end subroutine conservative_remap_3d

subroutine output_Tecplot_3d(nx, ny, nz, xGrid, yGrid, zGrid, u, v, w, T, filename)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: xGrid(nx), yGrid(ny), zGrid(nz)
    real(8), intent(in) :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), T(nx,ny,nz)
    character(len=*), intent(in) :: filename
    integer :: i, j, k
    real(4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5,V6,V7
    character(len=40) :: zoneName
    integer, parameter :: skip=1

    open(41, file=filename, access='stream', form='unformatted')
    zoneMarker = 299.0; eohMarker = 357.0
    
    ! Header Section
    write(41) "#!TDV101"
    write(41) 1
    title="3D_Data"
    call dumpstring(title)
    write(41) 7 ! Number of variables
    V1='X'; call dumpstring(V1)
    V2='Y'; call dumpstring(V2)
    V3='Z'; call dumpstring(V3)
    V4='U'; call dumpstring(V4)
    V5='V'; call dumpstring(V5)
    V6='W'; call dumpstring(V6)
    V7='T'; call dumpstring(V7)

    ! Zone Header
    write(41) zoneMarker
    zoneName="Zone 1"; call dumpstring(zoneName)
    write(41) -1
    write(41) 0 ! ZoneType = ORDERED
    write(41) 1 ! DataPacking = BLOCK
    write(41) 0
    write(41) 0
    write(41) 1+(nx-1)/skip
    write(41) 1+(ny-1)/skip
    write(41) 1+(nz-1)/skip
    write(41) 0
    write(41) eohMarker

    ! Zone Data
    write(41) zoneMarker
    write(41) 1; write(41) 1; write(41) 1 ! Data type: 1=float
    write(41) 1; write(41) 1; write(41) 1
    write(41) 1
    write(41) 0
    write(41) -1
    
    ! Write data in BLOCK format
    do k=1,nz,skip
        do j=1,ny,skip
            do i=1,nx,skip
                write(41) real(xGrid(i))
                write(41) real(yGrid(j))
                write(41) real(zGrid(k))
                write(41) real(u(i,j,k))
                write(41) real(v(i,j,k))
                write(41) real(w(i,j,k))
                write(41) real(T(i,j,k))
            end do
        end do
    enddo
    close(41)
    return

contains
    subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer :: stringLength
    integer :: ii
    integer :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
    end subroutine dumpstring
end subroutine output_Tecplot_3d