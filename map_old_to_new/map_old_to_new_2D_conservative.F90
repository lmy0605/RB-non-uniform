program map_old_to_new_2d_nonuniform
    implicit none

    !--------------------------- Variable Declarations ----------------------------
    character(len=100), parameter :: oldFile = "buoyancyCavity-1000.bin"
    character(len=100), parameter :: newFile = "u2u_g2l.bin"

    ! Grid dimensions
    integer, parameter :: nx_old = 257, ny_old = nx_old
    integer, parameter :: nx_new = 129, ny_new = nx_new

    ! Grid distribution types
    logical, parameter :: old_uniform = .true.  ! Set to .true. for uniform grid
    logical, parameter :: new_uniform = .true.  ! Set to .true. for uniform grid

    ! Constants for grid distribution
    real(8), parameter :: stretch_factor_old = 1.5d0
    real(8), parameter :: stretch_factor_new = 1.1d0

    ! Coordinate arrays for grids
    real(8), allocatable :: x_old(:), y_old(:)
    real(8), allocatable :: x_new(:), y_new(:)

    ! Data arrays for old and new grids
    real(8), allocatable :: u_old(:,:), v_old(:,:), T_old(:,:), rho_old(:,:)
    real(8), allocatable :: u_new(:,:), v_new(:,:), T_new(:,:), rho_new(:,:)

    integer :: ios, i, j

    integer, parameter :: read_unit = 10, write_unit = 20

    !--------------------------- Generate Mesh ----------------------------
    allocate(x_old(nx_old), y_old(ny_old))
    allocate(x_new(nx_new), y_new(ny_new))
    ! Compute nonuniform grid coordinates for old grid
    if (old_uniform .eqv. .true.) then
        do i = 1, nx_old
            x_old(i) = (dble(i) - 0.5d0) / dble(nx_old)
        end do
        do j = 1, ny_old
            y_old(j) = (dble(j) - 0.5d0) / dble(ny_old)
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
    end if

    ! Compute nonuniform grid coordinates for new grid
    if (new_uniform .eqv. .true.) then
        do i = 1, nx_new
            x_new(i) = (dble(i) - 0.5d0) / dble(nx_new)
        end do
        do j = 1, ny_new
            y_new(j) = (dble(j) - 0.5d0) / dble(ny_new)
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
    end if
    !--------------------------- Generate Mesh ----------------------------

    !--------------------------- Data Input ----------------------------
    allocate(u_old(nx_old,ny_old), v_old(nx_old,ny_old), T_old(nx_old,ny_old), rho_old(nx_old,ny_old))
    open(unit=read_unit, file=oldFile, form="unformatted", access="sequential", status='old', iostat=ios)    
    read(read_unit) ((u_old(i,j), i=1,nx_old), j=1,ny_old)
    read(read_unit) ((v_old(i,j), i=1,nx_old), j=1,ny_old)
    read(read_unit) ((T_old(i,j), i=1,nx_old), j=1,ny_old)
    ! read(read_unit) ((rho_old(i,j), i=1,nx_old), j=1,ny_old)
    !~!~rho_old = 1.0d0  ! Set density to a constant value for simplicity
    close(read_unit)

    !--------------------------- Check for Errors ----------------------------
    write(*,*) "old data loaded successfully!"
    if (old_uniform .eqv. .false.) then
        write(*,*) "old mesh is nonuniform"
        write(*,*) "stretch factor (old)=", real(stretch_factor_old)
    else
        write(*,*) "old mesh is uniform; NO stretch factor"
    end if
    write(*,*) "nx_old=", nx_old, " ny_old=", ny_old
    write(*,*) "u_old (center)=", real(u_old((nx_old-1)/2+1,(ny_old-1)/2+1))
    write(*,*) "v_old (center)=", real(v_old((nx_old-1)/2+1,(ny_old-1)/2+1))
    write(*,*) "T_old (center)=", real(T_old((nx_old-1)/2+1,(ny_old-1)/2+1))
    write(*,*) "rho_old (center)=", real(rho_old((nx_old-1)/2+1,(ny_old-1)/2+1))
    call output_Tecplot(nx_old, ny_old, x_old, y_old, u_old, v_old, T_old, "old.plt")
    !--------------------------- Check for Errors ----------------------------

    !--------------------------- Interpolation ----------------------------
    allocate(u_new(nx_new,ny_new), v_new(nx_new,ny_new), T_new(nx_new,ny_new), rho_new(nx_new,ny_new))
    call conservative_remap_2d( &
        x_old, y_old, u_old, v_old, T_old, rho_old, &
        x_new, y_new, u_new, v_new, T_new, rho_new, &
        nx_old, ny_old, nx_new, ny_new)

    !--------------------------- Check for Errors ----------------------------
    write(*,*) "Interpolation completed successfully!"
    if (new_uniform .eqv. .false.) then
        write(*,*) "new mesh is nonuniform"
        write(*,*) "stretch factor (new)=", real(stretch_factor_new)
    else
        write(*,*) "new mesh is uniform; NO stretch factor"
    end if
    write(*,*) "nx_new=", nx_new, " ny_new=", ny_new
    write(*,*) "u_new (center)=", real(u_new((nx_new-1)/2+1,(ny_new-1)/2+1))
    write(*,*) "v_new (center)=", real(v_new((nx_new-1)/2+1,(ny_new-1)/2+1))    
    write(*,*) "T_new (center)=", real(T_new((nx_new-1)/2+1,(ny_new-1)/2+1))
    write(*,*) "rho_new (center)=", real(rho_new((nx_new-1)/2+1,(ny_new-1)/2+1))
    call output_Tecplot(nx_new, ny_new, x_new, y_new, u_new, v_new, T_new, "new.plt")
    !--------------------------- Check for Errors ----------------------------

    !--------------------------- Data Output ----------------------------
    open(unit=write_unit, file=newFile, form="unformatted", access="sequential", status='replace', iostat=ios)
    write(write_unit) ((u_new(i,j), i=1,nx_new), j=1,ny_new)
    write(write_unit) ((v_new(i,j), i=1,nx_new), j=1,ny_new)
    write(write_unit) ((T_new(i,j), i=1,nx_new), j=1,ny_new)
    write(write_unit) ((rho_new(i,j), i=1,nx_new), j=1,ny_new)
    close(write_unit)

    write(*,*) "new data written successfully!"

    deallocate(x_old, y_old, x_new, y_new)
    deallocate(u_old, v_old, T_old, rho_old)
    deallocate(u_new, v_new, T_new, rho_new)

    print *, "2D mapping completed successfully!"
end program map_old_to_new_2d_nonuniform

subroutine conservative_remap_2d( &
    x_src, y_src, u_src, v_src, T_src, rho_src, &
    x_tgt, y_tgt, u_tgt, v_tgt, T_tgt, rho_tgt, &
    nx_src, ny_src, nx_tgt, ny_tgt)
    
    integer, intent(in) :: nx_src, ny_src, nx_tgt, ny_tgt
    real(8), intent(in) :: x_src(nx_src), y_src(ny_src)
    real(8), intent(in) :: u_src(nx_src, ny_src)
    real(8), intent(in) :: v_src(nx_src, ny_src)
    real(8), intent(in) :: T_src(nx_src, ny_src)
    real(8), intent(in) :: rho_src(nx_src, ny_src)
    real(8), intent(in) :: x_tgt(nx_tgt), y_tgt(ny_tgt)
    real(8), intent(out) :: u_tgt(nx_tgt, ny_tgt)
    real(8), intent(out) :: v_tgt(nx_tgt, ny_tgt)
    real(8), intent(out) :: T_tgt(nx_tgt, ny_tgt)
    real(8), intent(out) :: rho_tgt(nx_tgt, ny_tgt)

    ! --- Local variables ---
    integer :: i, j, ii, jj
    ! Dynamic search range indices
    integer :: i_start, i_end, j_start, j_end
    
    real(8) :: total_area
    real(8) :: u_integral, v_integral, T_integral, rho_integral
    real(8) :: overlap_area, x_overlap_center, y_overlap_center, u_recon, v_recon, T_recon, rho_recon

    ! --- Pre-calculated grid properties ---
    real(8), allocatable :: x_src_bounds(:,:), y_src_bounds(:,:)
    real(8), allocatable :: x_tgt_bounds(:,:), y_tgt_bounds(:,:)
    
    ! --- Gradients on source grid ---
    real(8), allocatable :: u_grad_x(:,:), u_grad_y(:,:)
    real(8), allocatable :: v_grad_x(:,:), v_grad_y(:,:)
    real(8), allocatable :: T_grad_x(:,:), T_grad_y(:,:)
    real(8), allocatable :: rho_grad_x(:,:), rho_grad_y(:,:)

    ! 1. Allocate space for gradients
    allocate(u_grad_x(nx_src, ny_src), u_grad_y(nx_src, ny_src))
    allocate(v_grad_x(nx_src, ny_src), v_grad_y(nx_src, ny_src))
    allocate(T_grad_x(nx_src, ny_src), T_grad_y(nx_src, ny_src))
    allocate(rho_grad_x(nx_src, ny_src), rho_grad_y(nx_src, ny_src))

    ! 2. Compute limited gradients for all source fields
    call compute_limited_gradients(u_src, x_src, y_src, nx_src, ny_src, u_grad_x, u_grad_y)
    call compute_limited_gradients(v_src, x_src, y_src, nx_src, ny_src, v_grad_x, v_grad_y)
    call compute_limited_gradients(T_src, x_src, y_src, nx_src, ny_src, T_grad_x, T_grad_y)
    call compute_limited_gradients(rho_src, x_src, y_src, nx_src, ny_src, rho_grad_x, rho_grad_y)

    ! 3. Pre-calculate cell boundaries to improve performance
    allocate(x_src_bounds(2, nx_src), y_src_bounds(2, ny_src))
    allocate(x_tgt_bounds(2, nx_tgt), y_tgt_bounds(2, ny_tgt))
    call calculate_cell_bounds(x_src, nx_src, x_src_bounds)
    call calculate_cell_bounds(y_src, ny_src, y_src_bounds)
    call calculate_cell_bounds(x_tgt, nx_tgt, x_tgt_bounds)
    call calculate_cell_bounds(y_tgt, ny_tgt, y_tgt_bounds)

    ! 4. Main loop over target grid
    do j = 1, ny_tgt
        do i = 1, nx_tgt
            ! Find the start and end indices of source cells that overlap with the current target cell.
            ! A buffer of -1 and +1 is added for safety, especially with non-uniform grids.
            i_start = max(1, bounded_bisect(x_src, nx_src, x_tgt_bounds(1, i)) - 1)
            i_end   = min(nx_src, bounded_bisect(x_src, nx_src, x_tgt_bounds(2, i)) + 1)
            j_start = max(1, bounded_bisect(y_src, ny_src, y_tgt_bounds(1, j)) - 1)
            j_end   = min(ny_src, bounded_bisect(y_src, ny_src, y_tgt_bounds(2, j)) + 1)
            
            ! Initialize integrals for the current target cell
            u_integral = 0.0d0; v_integral = 0.0d0
            T_integral = 0.0d0; rho_integral = 0.0d0
            total_area = 0.0d0

            ! Loop over the dynamically determined range
            do ii = i_start, i_end
                do jj = j_start, j_end
                    
                    ! Calculate the area and center of the overlapping region
                    call calculate_overlap_properties( &
                        x_src_bounds(:, ii), y_src_bounds(:, jj), &
                        x_tgt_bounds(:, i),  y_tgt_bounds(:, j), &
                        overlap_area, x_overlap_center, y_overlap_center)

                    if (overlap_area > 1.0d-14) then
                        ! Reconstruct value at the overlap center using source cell data and gradients
                        u_recon = u_src(ii,jj) + u_grad_x(ii,jj) * (x_overlap_center - x_src(ii)) &
                                               + u_grad_y(ii,jj) * (y_overlap_center - y_src(jj))
                        v_recon = v_src(ii,jj) + v_grad_x(ii,jj) * (x_overlap_center - x_src(ii)) &
                                               + v_grad_y(ii,jj) * (y_overlap_center - y_src(jj))
                        T_recon = T_src(ii,jj) + T_grad_x(ii,jj) * (x_overlap_center - x_src(ii)) &
                                               + T_grad_y(ii,jj) * (y_overlap_center - y_src(jj))
                        rho_recon = rho_src(ii,jj) + rho_grad_x(ii,jj) * (x_overlap_center - x_src(ii)) &
                                                 + rho_grad_y(ii,jj) * (y_overlap_center - y_src(jj))

                        ! Add the integrated contribution
                        u_integral = u_integral + u_recon * overlap_area
                        v_integral = v_integral + v_recon * overlap_area
                        T_integral = T_integral + T_recon * overlap_area
                        rho_integral = rho_integral + rho_recon * overlap_area
                        total_area = total_area + overlap_area
                    end if
                end do
            end do

            ! Normalize to get the final average value
            if (total_area > 1.0d-12) then
                u_tgt(i,j) = u_integral / total_area
                v_tgt(i,j) = v_integral / total_area
                T_tgt(i,j) = T_integral / total_area
                rho_tgt(i,j) = rho_integral / total_area
            else  
                ! Fallback to nearest neighbor (only if no overlap found, which is very unlikely with the fix)
                ii = bounded_bisect(x_src, nx_src, x_tgt(i))
                jj = bounded_bisect(y_src, ny_src, y_tgt(j))
                u_tgt(i,j) = u_src(ii,jj)
                v_tgt(i,j) = v_src(ii,jj)
                T_tgt(i,j) = T_src(ii,jj)
                rho_tgt(i,j) = rho_src(ii,jj)
            end if
        end do
    end do

    ! 5. Deallocate all temporary arrays
    deallocate(u_grad_x, u_grad_y, v_grad_x, v_grad_y, T_grad_x, T_grad_y, rho_grad_x, rho_grad_y)
    deallocate(x_src_bounds, y_src_bounds, x_tgt_bounds, y_tgt_bounds)

contains

    subroutine compute_limited_gradients(field, x, y, nx, ny, grad_x, grad_y)
        real(8), intent(in) :: field(nx, ny), x(nx), y(ny)
        integer, intent(in) :: nx, ny
        real(8), intent(out) :: grad_x(nx, ny), grad_y(nx, ny)
        
        integer :: i, j
        real(8) :: slope_L, slope_R, slope_B, slope_T, slope_cen

        ! Compute gradients in x-direction with van Leer limiter
        do j = 1, ny
            do i = 1, nx
                if (i > 1 .and. i < nx) then
                    slope_L = (field(i,j) - field(i-1,j)) / (x(i) - x(i-1))
                    slope_R = (field(i+1,j) - field(i,j)) / (x(i+1) - x(i))
                    slope_cen = (field(i+1,j) - field(i-1,j)) / (x(i+1) - x(i-1))
                    grad_x(i,j) = van_leer_limiter(slope_cen, slope_L, slope_R)
                else if (i == 1) then ! Forward difference
                    grad_x(i,j) = (field(2,j) - field(1,j)) / (x(2) - x(1))
                else ! Backward difference (i == nx)
                    grad_x(i,j) = (field(nx,j) - field(nx-1,j)) / (x(nx) - x(nx-1))
                end if
            end do
        end do
        
        ! Compute gradients in y-direction with van Leer limiter
        do i = 1, nx
            do j = 1, ny
                if (j > 1 .and. j < ny) then
                    slope_B = (field(i,j) - field(i,j-1)) / (y(j) - y(j-1))
                    slope_T = (field(i,j+1) - field(i,j)) / (y(j+1) - y(j))
                    slope_cen = (field(i,j+1) - field(i,j-1)) / (y(j+1) - y(j-1))
                    grad_y(i,j) = van_leer_limiter(slope_cen, slope_B, slope_T)
                else if (j == 1) then ! Forward difference
                    grad_y(i,j) = (field(i,2) - field(i,1)) / (y(2) - y(1))
                else ! Backward difference (j == ny)
                    grad_y(i,j) = (field(i,ny) - field(i,ny-1)) / (y(ny) - y(ny-1))
                end if
            end do
        end do
    end subroutine compute_limited_gradients

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
        real(8), intent(out) :: bounds(2, n) ! bounds(1,:) = left/bottom, bounds(2,:) = right/top
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

    subroutine calculate_overlap_properties(b_src_x, b_src_y, b_tgt_x, b_tgt_y, area, xc, yc)
        real(8), intent(in) :: b_src_x(2), b_src_y(2), b_tgt_x(2), b_tgt_y(2)
        real(8), intent(out) :: area, xc, yc
        real(8) :: overlap_dx, overlap_dy, x_min, x_max, y_min, y_max

        x_min = max(b_src_x(1), b_tgt_x(1))
        x_max = min(b_src_x(2), b_tgt_x(2))
        y_min = max(b_src_y(1), b_tgt_y(1))
        y_max = min(b_src_y(2), b_tgt_y(2))
        
        overlap_dx = x_max - x_min
        overlap_dy = y_max - y_min

        if (overlap_dx > 0.0d0 .and. overlap_dy > 0.0d0) then
            area = overlap_dx * overlap_dy
            xc = 0.5d0 * (x_min + x_max)
            yc = 0.5d0 * (y_min + y_max)
        else
            area = 0.0d0
            xc = 0.0d0
            yc = 0.0d0
        end if
    end subroutine calculate_overlap_properties

    function bounded_bisect(arr, n, x) result(idx)
        integer, intent(in) :: n
        real(8), intent(in) :: arr(n), x
        integer :: idx, lo, hi, mid
        
        if (x <= arr(1)) then
            idx = 1
            return
        else if (x >= arr(n)) then
            idx = n
            return
        end if

        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (x >= arr(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        idx = lo
    end function bounded_bisect

end subroutine conservative_remap_2d

subroutine output_Tecplot(nx, ny, xGrid, yGrid, u, v, T, filename)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: xGrid(nx), yGrid(ny)
    real(8), intent(in) :: u(nx, ny)
    real(8), intent(in) :: v(nx, ny)
    real(8), intent(in) :: T(nx, ny)
    character(len=*), intent(in) :: filename
    integer(kind=4) :: i, j, k
    REAL(kind=4) :: zoneMarker, eohMarker
    character(len=40) :: title
    character(len=40) :: V1,V2,V3,V4,V5
    integer(kind=4), parameter :: kmax=1
    character(len=40) :: zoneName

    open(41, file=filename, access='stream', form='unformatted')
    !---------------------------------------------
    zoneMarker= 299.0
    eohMarker = 357.0
    !I. HEAD SECTION--------------------------------------
    !c--Magic number, Version number
    write(41) "#!TDV101"

    !c--Integer value of 1
    write(41) 1

    Title="MyFirst"
    call dumpstring(title)
    !c-- Number of variables in this data file
    write(41) 5

    V1='X'
    call dumpstring(V1)
    V2='Y'
    call dumpstring(V2)
    V3='U'
    call dumpstring(V3)
    V4='V'
    call dumpstring(V4)
    V5='T'
    call dumpstring(V5)

    write(41) zoneMarker
    !--------Zone name.
    zoneName="ZONE 001"
    call dumpstring(zoneName)
    !---------Zone Color
    write(41) -1
    !---------ZoneType
    write(41) 0
    !---------DataPacking 0=Block, 1=Point
    write(41) 1
    !---------Specify Var Location. 0 = Do not specify, all data
    !---------is located at the nodes. 1 = Specify
    write(41) 0
    !---------Number of user defined face neighbor connections
    ! (value >= 0)
    write(41) 0
    !---------IMax,JMax,KMax
    write(41) nx
    write(41) ny
    write(41) kmax
    !-----------1=Auxiliary name/value pair to follow
    !-----------0=No more Auxiliar name/value pairs.
    write(41) 0
    write(41) eohMarker
    !----zone ------------------------------------------------------------
    write(41) zoneMarker
    !--------variable data format, 1=Float, 2=Double, 3=LongInt,4=ShortInt, 5=Byte, 6=Bit
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 1
    write(41) 0
    write(41) -1
    do k=1,kmax
        do j=1,ny
            do i=1,nx
                write(41) real(xGrid(i))
                write(41) real(yGrid(j))
                write(41) real(u(i,j))
                write(41) real(v(i,j))
                write(41) real(T(i,j))
            end do
        end do
    enddo
    close(41)
    return
    end subroutine output_Tecplot

    subroutine dumpstring(instring)
    implicit none
    character(len=40) instring
    integer(kind=4) :: stringLength
    integer(kind=4) :: ii
    integer(kind=4) :: I

    stringLength=LEN_TRIM(instring)
    do ii=1,stringLength
        I=ICHAR(instring(ii:ii))
        write(41) I
    end do
    write(41) 0

    return
    end subroutine dumpstring