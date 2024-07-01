module constant
    implicit none 

    real(8), parameter :: Pi = acos(-1.0d0)
    real(8), parameter :: g(0:2) = (/0.0d0, -1.0d0, 0.0d0/)
    real(8), parameter :: zero(0:2) = (/0.0d0, 0.0d0, 0.0d0/)

    !粒子のパラメータ
    integer, parameter :: N = 5000, Nstep = 500000
    real(8), parameter :: dt = 0.0005d0
    real(8), parameter :: m = 1.0d0, d1 = 1.0d0, d2 = 0.8d0
    real(8), parameter :: mu = 0.9d0, mu_cyliner = 0.9d0             !粒子, 円筒の密度
    real(8), parameter :: kn = 2.0d0 * 10 ** 5, kt = (2.0d0 / 7.0d0) * kn, e = 0.1d0
    real(8), parameter :: etan = -2 * log(e) * sqrt((m/2.0d0) * kn / ((acos(-1.0d0)) ** 2 + (log(e)) ** 2)), etat = etan


    !壁のパラメータ
    real(8), parameter :: Dwall = 50.0d0, W_wall(0:2) = (/0.0d0, 0.0d0, 0.1d0/)
    real(8), parameter :: etan_wall = -2 * log(e) * sqrt(m * kn / ((acos(-1.0d0)) ** 2 + (log(e)) ** 2)), etat_wall = etan_wall


    !linked-cell
    integer, parameter :: Nc = int(0.9d0*Dwall/d1)
    real(8), parameter :: cell_size = Dwall / Nc
    integer, parameter :: surplus_cell = 1
    integer, parameter :: Ncell = (Nc + 2*surplus_cell) ** 3

    integer, save :: icell_x = 0, icell_y = 0, icell_z = 0
    integer, save :: NextOf(1:N) = 0
    integer, save :: first(1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell) = 0
    integer, save :: last(1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell) = 0



    !グローバル変数
    real(8), save :: x_prev(0:2,1:N) = 0, v_prev(0:2,1:N) = 0
    real(8), save :: x(0:2,1:N) = 0, v(0:2,1:N) = 0
    real(8), save :: theta_prev(0:2,1:N) = 0, w_prev(0:2,1:N) = 0
    real(8), save :: theta(0:2,1:N) = 0, w(0:2,1:N) = 0
    real(8), save :: d(1:N) = 0, inertia(1:N) = 0
    real(8), save :: fn(0:2) = 0, ft(0:2) = 0
    real(8), save :: fall(0:2,1:N) = 0
    real(8), save :: deltat_prev(0:2,1:N*(N-1)/2) = 0
    real(8), save :: torque(0:2) = 0, torque_all(0:2,1:N) = 0
    real(8), save :: deltat_prev_wall(0:2,1:N) = 0.0d0

    real(8), save :: d_vec(0:2,1:N) = 0.0d0

end module constant

module mkfile
    use constant
    implicit none
contains
    subroutine makedirs(outdir)                                    !指定したディレクトリがなければ作る
        character(len=*), intent(in) :: outdir
        character(len=256) command
        write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
        ! write(*, *) trim(command)
        call system(command)
    end subroutine makedirs

    subroutine mkvtk(step, dirname)
        integer, intent(in) :: step
        character(len=*), intent(in) :: dirname
        integer i
        character(len=40) str_step
        character(len=40) str_N

        write(str_step, '(i7.7)') step
        write(str_N, '(i4)') N
        open(10,file = './'//trim(dirname)//'/data'//trim(str_step)//'.vtk')

        write(10,'(a)') '# vtk DataFile Version 5.1'
        write(10,'(a)') 'Data.vtk'
        write(10,'(a)') 'ASCII'
        write(10,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(10,'(a)') ' '
        write(10,'(a)') 'POINTS '//trim(str_N)//' float'
        do i = 1, N
            write(10,*) x(:,i)
        enddo
        write(10,'(a)') ' '
        write(10,'(a)') 'CELL_TYPES '//trim(str_N)
        do i = 1, N
            write(10,*) 1
        enddo
        write(10,'(a)') ' '
        write(10,'(a)') 'POINT_DATA '//trim(str_N)
        write(10,'(a)') 'SCALARS radius float'
        write(10,'(a)') 'LOOKUP_TABLE default'
        do i = 1, N 
            write(10,*) d(i)/2.0d0
        enddo

        close(10)
    end subroutine mkvtk

    subroutine mkvtk_sphere_wall(step, dirname)
        integer, intent(in) :: step
        character(len=*), intent(in) :: dirname
        integer i
        character(len=40) str_step
        character(len=40) str_N

        write(str_step, '(i7.7)') step
        write(str_N, '(i4)') 1
        open(10,file = './'//trim(dirname)//'/data_wall'//trim(str_step)//'.vtk')

        write(10,'(a)') '# vtk DataFile Version 5.1'
        write(10,'(a)') 'Data.vtk'
        write(10,'(a)') 'ASCII'
        write(10,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(10,'(a)') ' '
        write(10,'(a)') 'POINTS '//trim(str_N)//' float'
        write(10,*) 0, 0, 0
        write(10,'(a)') ' '
        write(10,'(a)') 'CELL_TYPES '//trim(str_N)
        write(10,*) 1
        write(10,'(a)') ' '
        write(10,'(a)') 'POINT_DATA '//trim(str_N)
        write(10,'(a)') 'SCALARS radius float'
        write(10,'(a)') 'LOOKUP_TABLE default'
        write(10,*) Dwall/2.0d0


        close(10)
    end subroutine mkvtk_sphere_wall



end module mkfile

module calculate 
    use constant
    implicit none 
contains 
    function vec_length(p, q) result(length)         !2点間の距離を計算する関数
        real(8), intent(in) :: p(:), q(:)
        real(8) length, sum
        integer i

        sum = 0
        do i = 1, size(p)
            sum = sum + (p(i) - q(i))**2
        enddo
        length = sqrt(sum)
    end function vec_length

    function inner_product(p, q) result(product)       !内積を計算するサブルーチン
        real(8), intent(in) :: p(:), q(:)
        real(8) product
        integer i 

        product = 0
        do i = 1, size(p)
            product = product + p(i) * q(i)
        enddo
    end function inner_product

    subroutine cross_product(p, q, cross)
        real(8), intent(in) :: p(0:2), q(0:2)
        real(8), intent(out) :: cross(0:2) 

        cross(0) = p(1)*q(2) - p(2)*q(1)
        cross(1) = p(2)*q(0) - p(0)*q(2)
        cross(2) = p(0)*q(1) - p(1)*q(0)
    end subroutine cross_product


    subroutine cal_direction_vector
        real(8) :: dtheta(0:2) = 0.0d0
        real(8) :: d_vec_tmp(0:2) = 0.0d0
        real(8) :: Rot_mat(0:2,0:2) = 0.0d0
        integer i

        do i = 1, N
            dtheta(:) = w(:,i) * dt
            Rot_mat(0,0) = cos(dtheta(2))*cos(dtheta(1))
            Rot_mat(0,1) = cos(dtheta(2))*sin(dtheta(1))*sin(dtheta(0))-cos(dtheta(0))*sin(dtheta(2))
            Rot_mat(0,2) = sin(dtheta(2))*sin(dtheta(0))+cos(dtheta(2))*cos(dtheta(0))*sin(dtheta(1))
            Rot_mat(1,0) = cos(dtheta(1))*sin(dtheta(2))
            Rot_mat(1,1) = cos(dtheta(2))*cos(dtheta(0))+sin(dtheta(2))*sin(dtheta(1))*sin(dtheta(0))
            Rot_mat(1,2) = cos(dtheta(0))*sin(dtheta(2))*sin(dtheta(1))-cos(dtheta(2))*sin(dtheta(0))
            Rot_mat(2,0) = -sin(dtheta(1))
            Rot_mat(2,1) = cos(dtheta(1))*sin(dtheta(0))
            Rot_mat(2,2) = cos(dtheta(1))*cos(dtheta(0))

            d_vec_tmp(0) = inner_product(Rot_mat(0,:),d_vec(:,i))
            d_vec_tmp(1) = inner_product(Rot_mat(1,:),d_vec(:,i))
            d_vec_tmp(2) = inner_product(Rot_mat(2,:),d_vec(:,i))

            d_vec(:,i) = d_vec_tmp(:)
            d_vec_tmp(:) = 0.0d0
            dtheta(:) = 0.0d0
        enddo

    end subroutine cal_direction_vector



    subroutine cal_cell_num(i)    !i番目の粒子がいるセル番号を計算するサブルーチン
        integer, intent(in) :: i
        icell_x = floor((x_prev(0,i) + Dwall/2.0d0) / cell_size) + 1 + surplus_cell    !x方向のセル番号
        icell_y = floor((x_prev(1,i) + Dwall/2.0d0) / cell_size) + 1 + surplus_cell    !y方向のセル番号
        icell_z = floor((x_prev(2,i) + Dwall/2.0d0) / cell_size) + 1 + surplus_cell    !z方向のセル番号
    end subroutine cal_cell_num


    subroutine mk_NextOf    !ある時刻でのNextOf配列を作成
        integer i
        integer last_prev

        first(:,:,:) = -1
        last(:,:,:) = -1
        NextOf(:) = -1

        do i = 1, N
            call cal_cell_num(i)

            last_prev = last(icell_x, icell_y, icell_z)    !icell番目のセル内の最後尾の粒子番号をとっておく
            last(icell_x, icell_y, icell_z) = i               !最後尾の粒子番号をiに更新
    
            if (last_prev == -1) then    !最後尾の粒子番号が-1(空)であれば最初の粒子番号をiとする
                first(icell_x, icell_y, icell_z) = i
            else 
                NextOf(last_prev) = i
            endif 
        enddo 

    end subroutine mk_NextOf

end module calculate

module contact 
    use constant
    use calculate
    implicit none
contains
    subroutine initial_volue     !粒子の初期位置を設定するサブルーチン
        integer :: i = 0
        real(8) :: rnd1, rnd2, rnd3
        integer(4), allocatable :: seed(:)
        integer(4) :: c, sz

        call random_seed( size = sz ) !コンパイラによってシードの配列サイズが異なるので，まずはシードの配列サイズを取得する．
        allocate( seed(sz) ) !得られたサイズでseedをallocateする．
        call random_seed( get = seed ) !デフォルトでセットされているシードを取得する．
        call system_clock( count = c ) !システム時刻を取得して，
        seed(1) = c                    !シードにその時刻を使用することで，毎回異なる乱数となるようにする．
        call random_seed( put = seed ) !組み直したシードを設定する．
        
        if(N > 0) then 
            do i = 1, N

                call random_number(rnd1)
                call random_number(rnd2)
                call random_number(rnd3)

                x(0,i) = (Dwall/2.0d0 - (d1/2.0d0)) * sqrt(rnd1) * sin(Pi*rnd2) * cos(2*Pi*rnd3)
                x(1,i) = (Dwall/2.0d0 - (d1/2.0d0)) * sqrt(rnd1) * sin(Pi*rnd2) * sin(2*Pi*rnd3)
                x(2,i) = (Dwall/2.0d0 - (d1/2.0d0)) * sqrt(rnd1) * cos(Pi*rnd2)
            enddo
        endif 

    end subroutine initial_volue


    subroutine particle_contact(i, j)    !iがjから受ける力
        integer, intent(in) :: i, j
        integer u
        real(8) :: norm(0:2) = 0, t(0:2) = 0
        real(8) :: v_tra(0:2) = 0, v_rot(0:2) = 0
        real(8) :: vn(0:2) = 0, vt(0:2) = 0, wn(0:2) = 0, wt(0:2) = 0
        real(8) :: deltan(0:2) = 0, deltat(0:2) = 0

        if ((d(i) + d(j))/2.0d0 - vec_length(x_prev(0:2,i),x_prev(0:2,j)) >= 0.0d0) then 
            u = 1
        else 
            u = 0
        endif 

        norm(0:2) = (x_prev(0:2,i) - x_prev(0:2,j))/vec_length(x_prev(0:2,i),x_prev(0:2,j))
        vn(0:2) = inner_product(v_prev(0:2,i) - v_prev(0:2,j),norm(0:2)) * norm(0:2)
        v_tra(0:2) = v_prev(0:2,i) - v_prev(0:2,j) - vn(0:2)
        wn(0:2) = inner_product(w_prev(0:2,i) - w_prev(0:2,j),norm(0:2)) * norm(0:2)
        wt(0:2) = w_prev(0:2,i) - w_prev(0:2,j) - wn(0:2)
        call cross_product(norm(0:2), (d(i)/2.0d0)*(w_prev(0:2,i) + w_prev(0:2,j)), v_rot(0:2))
        vt(0:2) = v_tra(0:2) + v_rot(0:2)

        if (vec_length(vt(0:2),zero(0:2)) == 0) then 
            t(0:2) = 0
        else 
            t(0:2) = vt(0:2) / vec_length(vt(0:2),zero(0:2))
        endif

        deltan(0:2) = -((d(i) + d(j))/2.0d0 - vec_length(x_prev(0:2,i),x_prev(0:2,j))) * norm(0:2)
        deltat(0:2) = u * (inner_product(deltat_prev(0:2,(j-1)*(j-2)/2 + i),t(0:2)) * t(0:2) + vt(0:2) * dt)

        fn(0:2) = u * (-kn * deltan(0:2) - etan * vn(0:2))
        ft(0:2) = u * (-kt * deltat(0:2) - etat * vt(0:2))

        if (vec_length(ft(0:2),zero(0:2)) > mu * vec_length(fn(0:2),zero(0:2))) then 
            ft(0:2) = mu * vec_length(fn(0:2),zero(0:2)) * ft(0:2) / vec_length(ft(0:2),zero(0:2))
            deltat_prev(0:2,(j-1)*(j-2)/2 + i) = inner_product(deltat_prev(0:2,(j-1)*(j-2)/2 + i), t(0:2)) * t(0:2)
        else 
            deltat_prev(0:2,(j-1)*(j-2)/2 + i) = deltat(0:2)
        endif

            
        call cross_product(-norm(0:2), ft(0:2), torque(0:2))    !iがjから受けるトルクの向き(大きさは後で決定する)

    end subroutine particle_contact

    subroutine sphere_wall_contact(i)
        integer, intent(in) :: i
        integer u_wall
        real(8) Lpw 
        real(8) norm(0:2), t(0:2), vn(0:2), vt(0:2)
        real(8) v_tra(0:2), v_rot(0:2)
        real(8) deltat_wall(0:2)
        real(8) R_r(0:2)
        real(8) :: b_wall = 0.0d0, ovrp = 0.0d0
        real(8) frot(0:2)

        Lpw = vec_length(x_prev(:,i),zero(:))

        if (Lpw >= (Dwall - d(i))/2.0d0) then
            u_wall = 1
        else 
            u_wall = 0
        endif 

        ovrp = u_wall*(Lpw - (Dwall - d(i))/2.0d0)
        b_wall = sqrt(d(i)*ovrp - ovrp**2)

        norm(0:2) = x_prev(:,i) / vec_length(x_prev(:,i), zero(:))
        vn(0:2) = inner_product(v_prev(:,i), norm(:)) * norm(0:2)
        v_tra(0:2) = v_prev(0:2,i) - vn(0:2)

        if (vec_length(vt(0:2),zero(0:2)) == 0) then 
            t(0:2) = 0
        else 
            t(0:2) = vt(0:2) / vec_length(vt(0:2),zero(0:2))
        endif

        call cross_product(-norm(0:2), (d(i)/2.0d0) * w_prev(0:2,i) - (Dwall/2.0d0) * W_wall(0:2), v_rot(0:2))
        vt(:) = v_tra(:) + v_rot(:)
        deltat_wall(0:2) = u_wall * (inner_product(deltat_prev_wall(0:2,i),t(0:2)) * t(0:2) + vt(0:2) * dt)

        fn(:) = u_wall * (-kn * (Lpw - (Dwall - d(i))/2.0d0) * norm(:) - etan_wall * vn(:))
        ft(:) = u_wall * (-kt * deltat_wall(:) - etat_wall * vt(:))
        frot(:) = u_wall * (- etat_wall * v_rot(:))

        if (vec_length(ft(0:2),zero(0:2)) > mu * vec_length(fn(0:2),zero(0:2))) then 
            ft(0:2) = mu * vec_length(fn(0:2),zero(0:2)) * ft(0:2) / vec_length(ft(0:2),zero(0:2))
            deltat_prev_wall(0:2,i) = inner_product(deltat_prev_wall(0:2,i), t(0:2)) * t(0:2)
        else 
            deltat_prev_wall(0:2,i) = deltat_wall(0:2)
        endif

        R_r(:) = frot(:) * 2.0d0 * b_wall**2 / (3.0d0*Dwall)
            
        call cross_product(norm(0:2), ft(0:2), torque(0:2)) 
        torque(:) = torque(:) + R_r(:)
    end subroutine sphere_wall_contact
end module contact 

program dem_3d
    use constant
    use mkfile
    use calculate
    use contact
    implicit none 

    integer i, j, step
    integer ix, iy, iz
    integer :: progress = 0
    character(len=40) :: dirname = 'vtk_data'

    call makedirs(dirname)

    d(1:N/2) = d1
    d(N/2 + 1:N) = d2
    inertia(1:N/2) = 0.40d0 * m * ((d1 / 2.0d0) ** 2)
    inertia(N/2 + 1:N) = 0.40d0 * m * ((d2 / 2.0d0) ** 2)


    open(20, file = 'xy1')
    open(30, file = 'xy2')


    ! x(:,1) = (/1.0d0, 0.0d0, 0.0d0/)
    ! ! x(:,2) = (/-1.0d0, 0.0d0, 0.0d0/)
    ! v(:,1) = (/-1.0d0, 0.0d0, 0.0d0/)
    ! ! v(:,2) = (/1.0d0, 0.0d0, 0.0d0/)
    ! w(:,1) = (/0.0d0, 0.0d0, -1.0d0/)
    call initial_volue

    do i = 1, N
        d_vec(:,i) = (/0.0, 1.0, 0.0/)
    enddo

    call mkvtk(0, dirname)
    call mkvtk_sphere_wall(0, dirname)

    do step = 1, Nstep
        v_prev(:,:) = v(:,:)
        x_prev(:,:) = x(:,:)
        w_prev(:,:) = w(:,:)
        theta_prev(:,:) = theta(:,:)

        call mk_NextOf

        ! do i = 1, N - 1
        !     do j = i + 1, N
        !             !粒子間の接触
        !             call particle_contact(i, j)
        !             fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
        !             fall(0:2,j) = fall(0:2,j) - fn(0:2) - ft(0:2)
        !             torque_all(0:2,i) = torque_all(0:2,i) + d(i) * torque(0:2) / 2.0d0   !ここでトルクの大きさを決定
        !             torque_all(0:2,j) = torque_all(0:2,j) + d(j) * torque(0:2) / 2.0d0
        !             fn(:) = 0
        !             ft(:) = 0
        !             torque(:) = 0
        !     enddo 
        ! enddo

        do i = 1, N-1
            call cal_cell_num(i)
            do iz = -1, 1, 1
                do iy = -1, 1, 1
                    do ix = -1, 1, 1
                        j = first(icell_x+ix, icell_y+iy, icell_z+iz)

                        do 
                            if (j == -1) then 
                                exit
                            elseif (i < j) then
                                !粒子間の接触
                                call particle_contact(i, j)
                                fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
                                fall(0:2,j) = fall(0:2,j) - fn(0:2) - ft(0:2)
                                torque_all(0:2,i) = torque_all(0:2,i) + d(i) * torque(0:2) / 2.0d0   !ここでトルクの大きさを決定
                                torque_all(0:2,j) = torque_all(0:2,j) + d(j) * torque(0:2) / 2.0d0
                                fn(:) = 0
                                ft(:) = 0
                                torque(:) = 0
                            endif
                            j = nextOf(j)
                        enddo
                    enddo
                enddo
            enddo
        enddo




        do i = 1, N
            !壁との接触
            call sphere_wall_contact(i)
            fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
            torque_all(:,i) = torque_all(:,i) + d(i) * torque(0:2) / 2.0d0
            fn(:) = 0
            ft(:) = 0
            torque(:) = 0
            
            v(:,i) = v_prev(:,i) + (fall(:,i)/m + g(:)) * dt
            x(:,i) = x_prev(:,i) + v(:,i) * dt
            w(:,i) = w_prev(:,i) + torque_all(:,i) * dt / inertia(i)
            theta(:,i) = theta_prev(:,i) + w(:,i) * dt

            fall(:,i) = 0.0d0
            torque_all(:,i) = 0.0d0
        enddo 

        call cal_direction_vector

        if (mod(step, 500) == 0) then
            call mkvtk(step, dirname)
            call mkvtk_sphere_wall(step, dirname)
        endif

        write(20,*) step, w(2,1), w(2,2)
        write(30,*) x(0,2), x(1,2)

        if(step == (Nstep / 10) * (progress + 1)) then 
            progress = progress + 1
            print *, '進捗',progress * 10,'%'
        endif 
    enddo


    close(20)
    close(30)
end program dem_3d