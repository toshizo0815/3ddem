module constant
    implicit none 

    real(8), parameter :: Pi = acos(-1.0d0)
    ! real(8), parameter :: g(0:2) = (/0.0d0, -1.0d0, 0.0d0/)
    real(8), parameter :: zero(0:2) = (/0.0d0, 0.0d0, 0.0d0/)

    !もとのパラメータ
    real(8), parameter :: rho_particle = 3.6d0    !アルミナ
    ! real(8), parameter :: rho_particle = 2.5d0    !ガラス

    real(8), parameter :: rho_cylinder = 1.2d0   !円筒の密度

    !real(8), parameter :: d_particle = 1.0d0
    real(8), parameter :: d_particle = 2.0d0
    !real(8), parameter :: d_particle = 6.0d0

    !粒子のパラメータ
    integer, parameter :: N = 2000, Nstep = 500000, step_fixed = 100000
    real(8), parameter :: dt = 0.0005d0
    real(8), parameter :: m = 1.0d0, d1 = 1.0d0, d2 = 0.8d0
    real(8), parameter :: mu = 0.9d0, mu_cyliner = 0.9d0             !粒子, 円筒の密度
    real(8), parameter :: kn = 2.0d0 * 10 ** 5, kt = (2.0d0 / 7.0d0) * kn, e = 0.1d0
    real(8), parameter :: etan = -2 * log(e) * sqrt((m/2.0d0) * kn / ((acos(-1.0d0)) ** 2 + (log(e)) ** 2)), etat = etan

    !円筒のパラメータ
    real(8), parameter :: Dwall = 60.0d0 / d_particle, deltad = 3.0d0 / d_particle
    real(8), parameter :: H_bottom = 20.0d0 / d_particle    !原点を中心に左右にH_bottomだけとる
    real(8), parameter :: deltah = 1.0d0 / d_particle
    real(8), parameter :: m_cylinder = 3*(rho_cylinder/rho_particle)*((H_bottom+deltah)*(Dwall+2*deltad)**2-H_bottom*Dwall**2)
    real(8), parameter :: inertia_cylinder = (Pi*rho_cylinder/16.0d0)*((H_bottom+deltah)*(Dwall+2*deltad)**4-H_bottom*Dwall**4)
    real(8), parameter :: etan_wall = -2 * log(e) * sqrt(m * kn / ((acos(-1.0d0)) ** 2 + (log(e)) ** 2)), etat_wall = etan_wall

    real(8), parameter :: theta_slope_degree = 10.0d0
    real(8), parameter :: theta_slope = theta_slope_degree * Pi / 180.0d0
    real(8), parameter :: etan_slope = -2 * log(e) * sqrt(m_cylinder * kn / (Pi ** 2 + (log(e)) ** 2)), etat_slope = etan_slope

    real(8), parameter :: g(0:2) = (/sin(theta_slope), -cos(theta_slope), 0.0d0/)

    !linked-cell
    integer, parameter :: Nc = int(0.9d0*Dwall/d1)
    real(8), parameter :: cell_size = Dwall / Nc
    integer, parameter :: surplus_cell = 1
    integer, parameter :: Ncell = (Nc + 2*surplus_cell) ** 3

    integer, save :: icell_x = 0, icell_y = 0, icell_z = 0
    integer, save :: NextOf(1:N+1) = 0
    integer, save :: first(1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell) = 0
    integer, save :: last(1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell, 1:Nc+2*surplus_cell) = 0



    !グローバル変数
    real(8), save :: x_prev(0:2,1:N+1) = 0.0d0, v_prev(0:2,1:N+1) = 0.0d0
    real(8), save :: x(0:2,1:N+1) = 0.0d0, v(0:2,1:N+1) = 0.0d0
    real(8), save :: theta_prev(0:2,1:N+1) = 0.0d0, w_prev(0:2,1:N+1) = 0.0d0
    real(8), save :: theta(0:2,1:N+1) = 0.0d0, w(0:2,1:N+1) = 0.0d0

    real(8), save :: x_prev_wall(0:2) = 0.0d0, v_prev_wall(0:2) = 0.0d0
    real(8), save :: x_wall(0:2) = 0.0d0, v_wall(0:2) = 0.0d0
    real(8), save :: theta_prev_wall(0:2) = 0.0d0, w_prev_wall(0:2) = 0.0d0
    real(8), save :: theta_wall(0:2) = 0.0d0, w_wall(0:2) = 0.0d0

    real(8), save :: x_slope(0:2) = (/0.0d0, -(Dwall/2.0d0 + deltad), 0.0d0/), v_slope(0:2) = 0.0d0
    real(8), save :: x_slope_prev(0:2) = 0.0d0, v_slope_prev(0:2) = 0.0d0

    real(8), save :: alpha(0:2) = 0.0d0

    real(8), save :: d(1:N+1) = 0.0d0, inertia(1:N+1) = 0.0d0
    real(8), save :: fn(0:2) = 0.0d0, ft(0:2) = 0.0d0, ft_wall(0:2) = 0.0d0
    real(8), save :: fall(0:2,1:N+1) = 0.0d0, fall_wall(0:2) = 0.0d0
    real(8), save :: deltat_prev(0:2,1:N+1,1:N+1) = 0.0d0
    real(8), save :: torque(0:2) = 0.0d0, torque_all(0:2,1:N+1) = 0.0d0
    real(8), save :: torque_wall(0:2) = 0.0d0, torque_all_wall(0:2) = 0.0d0
    real(8), save :: deltat_prev_wall(0:2,1:N+1) = 0.0d0
    real(8), save :: deltat_prev_bottom_wall(0:2,1:N+1,1:2) = 0.0d0    !1は左の壁，2は右の壁
    real(8), save :: deltat_prev_slope(0:2) = 0.0d0

    real(8), save :: d_vec(0:2,1:N+1) = 0.0d0
    real(8), save :: cylinder_vec(0:2) = 0.0d0

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

        !点の位置情報
        write(10,'(a)') 'POINTS '//trim(str_N)//' float'
        do i = 1, N
            write(10,'(f8.2,f8.2,f8.2)') x(:,i)
        enddo
        write(10,'(a)') ' '

        !半径の情報
        write(10,'(a)') 'POINT_DATA '//trim(str_N)
        write(10,'(a)') 'SCALARS radius float'
        write(10,'(a)') 'LOOKUP_TABLE default'
        do i = 1, N 
            write(10,*) d(i)/2.0d0
        enddo
        write(10,'(a)') ' '

        !方向指定ベクトルの向き
        write(10,'(a)') 'VECTORS direction float'
        do i = 1, N 
            write(10,'(f6.2,f6.2,f6.2)') d_vec(:,i) * d(i)/2.0d0
        enddo

        close(10)
    end subroutine mkvtk


    subroutine mkvtk_cylinder(step, dirname)
        integer, intent(in) :: step
        character(len=*), intent(in) :: dirname
        character(len=40) str_step
        character(len=40) str_N

        write(str_step, '(i7.7)') step
        write(str_N, '(i4)') 2
        open(10,file = './'//trim(dirname)//'/data_cylinder'//trim(str_step)//'.vtk')

        write(10,'(a)') '# vtk DataFile Version 5.1'
        write(10,'(a)') 'Data.vtk'
        write(10,'(a)') 'ASCII'
        write(10,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(10,'(a)') ' '
        write(10,'(a)') 'POINTS '//trim(str_N)//' float'
        write(10,*) 0, 0, 0
        write(10,*) 0, 0, 0
        write(10,'(a)') ' '
        write(10,'(a)') 'CELL_TYPES '//trim(str_N)
        write(10,*) 1
        write(10,*) 1
        write(10,'(a)') ' '
        write(10,'(a)') 'POINT_DATA '//trim(str_N)
        write(10,'(a)') 'SCALARS radius float'
        write(10,'(a)') 'LOOKUP_TABLE default'
        write(10,*) Dwall/2.0d0
        write(10,*) Dwall/2.0d0 + deltad
        write(10,'(a)') ' '

        !方向指定ベクトルの向き
        write(10,'(a)') 'VECTORS direction float'
        write(10,'(f6.2,f6.2,f6.2)') (Dwall/2.0d0) * cylinder_vec(0:2)
        write(10,'(f6.2,f6.2,f6.2)') (Dwall/2.0d0 + deltad) * cylinder_vec(0:2)


        close(10)
    end subroutine mkvtk_cylinder



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

    function inner_product(p, q) result(product)       !内積を計算する関数
        real(8), intent(in) :: p(:), q(:)
        real(8) product
        integer i 

        product = 0
        do i = 1, size(p)
            product = product + p(i) * q(i)
        enddo
    end function inner_product

    subroutine cross_product(p, q, cross)            !外積を計算するサブルーチン
        real(8), intent(in) :: p(0:2), q(0:2)
        real(8), intent(out) :: cross(0:2) 

        cross(0) = p(1)*q(2) - p(2)*q(1)
        cross(1) = p(2)*q(0) - p(0)*q(2)
        cross(2) = p(0)*q(1) - p(1)*q(0)
    end subroutine cross_product


    subroutine cal_direction_vector                  !粒子の向きを表すベクトルを計算するサブルーチン
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
        deltat(0:2) = u * (inner_product(deltat_prev(0:2,i,j),t(0:2)) * t(0:2) + vt(0:2) * dt)

        fn(0:2) = u * (-kn * deltan(0:2) - etan * vn(0:2))
        ft(0:2) = u * (-kt * deltat(0:2) - etat * vt(0:2))

        if (vec_length(ft(0:2),zero(0:2)) > mu * vec_length(fn(0:2),zero(0:2))) then 
            ft(0:2) = mu * vec_length(fn(0:2),zero(0:2)) * ft(0:2) / vec_length(ft(0:2),zero(0:2))
            deltat_prev(0:2,i,j) = inner_product(deltat_prev(0:2,i,j), t(0:2)) * t(0:2)
        else 
            deltat_prev(0:2,i,j) = deltat(0:2)
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

        call cross_product(-norm(0:2), (d(i)/2.0d0) * w_prev(0:2,i) - (Dwall/2.0d0) * w_prev_wall(0:2), v_rot(0:2))
        vt(:) = v_tra(:) + v_rot(:)

        if (vec_length(vt(0:2),zero(0:2)) == 0) then 
            t(0:2) = 0
        else 
            t(0:2) = vt(0:2) / vec_length(vt(0:2),zero(0:2))
        endif

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

    subroutine cylinder_wall_contact(i)
        integer, intent(in) :: i
        integer u_wall
        real(8) Lpw 
        real(8) norm(0:2), t(0:2), vn(0:2), vt(0:2)
        real(8) t_wall(0:2)
        real(8) v_tra(0:2), v_rot(0:2)
        real(8) deltat_wall(0:2)
        real(8) R_r(0:2)
        real(8) :: b_wall = 0.0d0, ovrp = 0.0d0
        real(8) frot(0:2)
        real(8) x_project(0:2)

        x_project(:) = (/x_prev(0,i), x_prev(1,i), 0.0d0/)

        Lpw = vec_length(x_project(:), zero(:))

        if (Lpw >= (Dwall - d(i))/2.0d0) then
            u_wall = 1
        else 
            u_wall = 0
        endif 

        ovrp = u_wall*(Lpw - (Dwall - d(i))/2.0d0)
        ! b_wall = sqrt(d(i)*ovrp - ovrp**2)

        norm(0:2) = x_project(:) / vec_length(x_project(:), zero(:))
        t_wall(0:2) = (/-norm(1), norm(0), norm(2)/)
        vn(0:2) = inner_product(v_prev(:,i) - v_prev_wall(:), norm(:)) * norm(0:2)
        v_tra(0:2) = v_prev(0:2,i) - v_prev_wall(0:2) - vn(0:2)

        call cross_product(-norm(0:2), (d(i)/2.0d0) * w_prev(0:2,i) - (Dwall/2.0d0) * w_prev_wall(0:2), v_rot(0:2))
        vt(:) = v_tra(:) + v_rot(:)

        if (vec_length(vt(0:2),zero(0:2)) == 0) then 
            t(0:2) = 0
        else 
            t(0:2) = vt(0:2) / vec_length(vt(0:2),zero(0:2))
        endif

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

        ft_wall(0:2) = inner_product(ft(:),t_wall(:))*t_wall(0:2)

        ! R_r(:) = frot(:) * 2.0d0 * b_wall**2 / (3.0d0*Dwall)
            
        call cross_product(norm(0:2), ft(0:2), torque(0:2)) 
        call cross_product(norm(0:2), ft_wall(0:2), torque_wall(0:2))
        torque(:) = torque(:) !+ R_r(:)
    end subroutine cylinder_wall_contact

    subroutine cylinder_bottom_wall_contact(i,LorR)       !LorRは1(left)か2(right)
        integer, intent(in) :: i, LorR
        integer u_wall
        real(8) Lpw 
        real(8) x_projection(0:2)
        real(8) norm(0:2), t(0:2), vn(0:2), vt(0:2)
        real(8) v_tra(0:2), v_rot(0:2)
        real(8) v_rot_wall(0:2)
        real(8) deltat_wall(0:2)
        real(8) R_r(0:2)
        real(8) :: b_wall = 0.0d0, ovrp = 0.0d0
        real(8) frot(0:2)

        x_projection(0:2) = (/x_prev(0,i),x_prev(1,i),0.0d0/)

        Lpw = (-1.0d0)**(LorR+1) * x_prev(2,i) + H_bottom

        if (Lpw < d(i)/2.0d0) then
            u_wall = 1
        else 
            u_wall = 0
        endif 

        ovrp = u_wall*(d(i)/2.0d0 - Lpw)
        ! b_wall = sqrt(d(i)*ovrp - ovrp**2)

        norm(0:2) = (/0.0d0, 0.0d0, (-1.0d0)**(LorR)/)
        vn(0:2) = inner_product(v_prev(:,i), norm(:)) * norm(0:2)
        v_tra(0:2) = v_prev(0:2,i) - vn(0:2)

        call cross_product(-norm(0:2), (d(i)/2.0d0) * w_prev(0:2,i), v_rot(0:2))
        call cross_product(w_prev_wall(:), x_projection(:), v_rot_wall(:))
        vt(:) = v_tra(:) + v_rot(:) - v_rot_wall(:)

        if (vec_length(vt(0:2),zero(0:2)) == 0) then 
            t(0:2) = 0
        else 
            t(0:2) = vt(0:2) / vec_length(vt(0:2),zero(0:2))
        endif

        deltat_wall(0:2) = u_wall * (inner_product(deltat_prev_bottom_wall(0:2,i,LorR),t(0:2)) * t(0:2) + vt(0:2) * dt)

        fn(:) = u_wall * (-kn * ovrp * norm(:) - etan_wall * vn(:))
        ft(:) = u_wall * (-kt * deltat_wall(:) - etat_wall * vt(:))
        frot(:) = u_wall * (- etat_wall * v_rot(:))

        if (vec_length(ft(0:2),zero(0:2)) > mu * vec_length(fn(0:2),zero(0:2))) then 
            ft(0:2) = mu * vec_length(fn(0:2),zero(0:2)) * ft(0:2) / vec_length(ft(0:2),zero(0:2))
            deltat_prev_bottom_wall(0:2,i,LorR) = inner_product(deltat_prev_bottom_wall(0:2,i,LorR), t(0:2)) * t(0:2)
        else 
            deltat_prev_bottom_wall(0:2,i,LorR) = deltat_wall(0:2)
        endif

        ! R_r(:) = frot(:) * 2.0d0 * b_wall**2 / (3.0d0*Dwall)
            
        call cross_product(norm(0:2), ft(0:2), torque(0:2)) 
        torque(:) = torque(:) !+ R_r(:)
    end subroutine cylinder_bottom_wall_contact

    subroutine slope_contact
        integer u_wall
        real(8) norm(0:2), t(0:2), vn(0:2), vt(0:2)
        real(8) v_tra(0:2), v_rot(0:2)
        real(8) deltat_slope(0:2)
        real(8) R_r(0:2)
        real(8) :: b_wall = 0.0d0, ovrp = 0.0d0
        real(8) frot(0:2)

        if ((Dwall/2.0d0 + deltad) + x_slope_prev(1) > 0) then
            u_wall = 1
        else 
            u_wall = 0
        endif 

        ovrp = u_wall*((Dwall/2.0d0 + deltad) + x_slope_prev(1))
        ! b_wall = sqrt(Dwall*ovrp - ovrp**2)

        norm(0:2) = (/0.0d0, -1.0d0, 0.0d0/)
        t(0:2) = (/1.0d0, 0.0d0, 0.0d0/)

        vn(0:2) = inner_product(-v_slope_prev(:), norm(:)) * norm(0:2)
        v_tra(0:2) = inner_product(-v_slope_prev(:), t(:)) * t(0:2)

        call cross_product(-norm(0:2), (Dwall/2.0d0 + deltad) * w_prev_wall(0:2), v_rot(0:2))
        vt(:) = v_tra(:) + v_rot(:)


        deltat_slope(0:2) = u_wall * (inner_product(deltat_prev_slope(0:2),t(0:2)) * t(0:2) + vt(0:2) * dt)

        fn(:) = u_wall * (-kn * ovrp * norm(:) - etan_slope * vn(:))
        ft(:) = u_wall * (-kt * deltat_slope(:) - etat_slope * vt(:))
        frot(:) = u_wall * (- etat_wall * v_rot(:))

        if (vec_length(ft(0:2),zero(0:2)) > mu * vec_length(fn(0:2),zero(0:2))) then 
            ft(0:2) = mu * vec_length(fn(0:2),zero(0:2)) * ft(0:2) / vec_length(ft(0:2),zero(0:2))
            deltat_prev_slope(0:2) = inner_product(deltat_prev_slope(0:2), t(0:2)) * t(0:2)
        else 
            deltat_prev_slope(0:2) = deltat_slope(0:2)
        endif

        ! R_r(:) = frot(:) * 2.0d0 * b_wall**2 / (3.0d0*Dwall)
            
        call cross_product(norm(0:2), ft(0:2), torque(0:2)) 

    end subroutine slope_contact


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
    character(len=40) :: dirname = 'vtk_data_mb'

    call makedirs(dirname)

    d(1:N/2) = d1
    d(N/2 + 1:N) = d2
    inertia(1:N/2) = 0.40d0 * m * ((d1 / 2.0d0) ** 2)
    inertia(N/2 + 1:N) = 0.40d0 * m * ((d2 / 2.0d0) ** 2)


    open(20, file = 'xt_cylinder')


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
    call mkvtk_cylinder(0, dirname)

    do step = 1, Nstep
        v_prev(:,:) = v(:,:)
        x_prev(:,:) = x(:,:)
        w_prev(:,:) = w(:,:)
        theta_prev(:,:) = theta(:,:)
        w_prev_wall(:) = w_wall(:)
        theta_prev_wall(:) = theta_wall(:)
        v_slope_prev(:) = v_slope(:)
        x_slope_prev(:) = x_slope(:)

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
            !球壁との接触
            !call sphere_wall_contact(i)
            ! fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
            ! torque_all(:,i) = torque_all(:,i) + d(i) * torque(0:2) / 2.0d0
            ! fn(:) = 0
            ! ft(:) = 0
            ! torque(:) = 0

            !円筒側面との接触
            call cylinder_wall_contact(i)
            fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
            fall_wall(0:2) = fall_wall(0:2) - fn(0:2) - ft_wall(0:2)
            torque_all(:,i) = torque_all(:,i) + d(i) * torque(0:2) / 2.0d0
            torque_all_wall(:) = torque_all_wall(:) - Dwall * torque_wall(0:2) / 2.0d0
            fn(:) = 0.0d0
            ft(:) = 0.0d0
            torque(:) = 0.0d0
            ft_wall(:) = 0.0d0
            torque_wall(:) = 0.0d0

            !円筒底面との接触
            call cylinder_bottom_wall_contact(i,1)
            fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
            torque_all(:,i) = torque_all(:,i) + d(i) * torque(0:2) / 2.0d0
            fn(:) = 0
            ft(:) = 0
            torque(:) = 0

            call cylinder_bottom_wall_contact(i,2)
            fall(0:2,i) = fall(0:2,i) + fn(0:2) + ft(0:2)
            torque_all(:,i) = torque_all(:,i) + d(i) * torque(0:2) / 2.0d0
            fn(:) = 0
            ft(:) = 0
            torque(:) = 0
        enddo 

        if(step > step_fixed) then      !step = step_fixedとなるまで円筒を固定する
            call slope_contact
            fall_wall(:) = fall_wall(:) + fn(:) + ft(:)
            torque_all_wall(:) = torque_all_wall(:) + (Dwall/2.0d0 + deltad)*torque(:)
            fn(:) = 0
            ft(:) = 0
            torque(:) = 0

            alpha(:) = fall_wall(:)/m_cylinder + g(:)   !座標系の加速度
        endif

        do i = 1, N
            v(:,i) = v_prev(:,i) + (fall(:,i)/m + g(:) - alpha(:)) * dt
            x(:,i) = x_prev(:,i) + v(:,i) * dt
            w(:,i) = w_prev(:,i) + torque_all(:,i) * dt / inertia(i)
            theta(:,i) = theta_prev(:,i) + w(:,i) * dt

            fall(:,i) = 0.0d0
            torque_all(:,i) = 0.0d0
        enddo

        if(step > step_fixed) then       !step = step_fixedとなるまで円筒を固定する
            v_slope(:) = v_slope_prev(:) - alpha(:) * dt
            x_slope(:) = x_slope_prev(:) + v_slope(:) * dt
            w_wall(:) = w_prev_wall(:) + torque_all_wall(:) * dt / inertia_cylinder
            theta_wall(:) = theta_prev_wall(:) + w_wall(:) * dt
        endif

        cylinder_vec(0:2) = (/cos(theta_wall(2)), sin(theta_wall(2)), 0.0d0/)


        fall_wall(:) = 0.0d0
        torque_all_wall(:) = 0.0d0

        call cal_direction_vector

        if (mod(step, 500) == 0) then
            call mkvtk(step, dirname)
            call mkvtk_cylinder(step, dirname)
        endif

        write(20,*) step, -x_slope(:)

        if(step == (Nstep / 10) * (progress + 1)) then 
            progress = progress + 1
            print *, '進捗',progress * 10,'%'
        endif 
    enddo

    ! open(30,file='xy_cylinder_only')
    ! do step = 0, Nstep
    !     write(30,*) step, (0.5d0*((step*dt)**2)*m_cylinder*(sin(theta_slope))*(Dwall/2.0d0+deltad)**2)&
    !     /(inertia_cylinder+m_cylinder*(Dwall/2.0d0+deltad)**2)
    ! enddo
    ! close(30)

    close(20)
end program dem_3d