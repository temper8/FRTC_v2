module driven_current_module
    !! Driven Current Module
    use kind_module
    implicit none
    type DrivenCurrentResult
        real(wp) :: cup, cp
        !!
        real(wp) :: cum, cm
        !!
        real(wp) :: cup0, cp0
        !!
        real(wp) :: cum0, cm0
        !!

    contains
        procedure :: print => driven_current_result_print
        procedure :: save => driven_current_result_save

    end type DrivenCurrentResult

contains

    subroutine driven_current_result_print(this, time)
        class(DrivenCurrentResult), intent(in) :: this
        real(wp), intent(in)  :: time
        print *, '------- driven current ---------'
        print *, 'time=',time
        print *, 'cup=',this%cup,' cp=',this%cp
        print *, 'cum=',this%cum,' cm=',this%cm
        print *, 'cup0=',this%cup0,' cp0=',this%cp0
        print *, 'cum0=',this%cum0,' cm0=',this%cm0
        print *, 'sigma driven current, MA=', this%cp0 + this%cm0
        print *, 'driven current, MA=', this%cup + this%cum
        print *, '--------------------------------'
    end subroutine driven_current_result_print   

    subroutine driven_current_result_save(this, time)
        class(DrivenCurrentResult), intent(in) :: this
        real(wp), intent(in)  :: time        
        logical, save :: first_time=.TRUE.
        character(132) FNAME
        integer :: io
        FNAME = "lhcd/dc_result.dat"

		if (first_time) then
            open(newunit=io, file=FNAME, status="replace", action="write")
			write(io,'(100A22)'), "Time", 'cup', 'cp', 'cum', 'cm', 'cup0', 'cp0', 'cum0', 'cm0'
            close(io)
			first_time = .FALSE.
        else
            open(newunit=io, file=FNAME, position="append")
            write(io,'(100f22.14)') time,  this%cup, this%cp, this%cum, this%cm, this%cup0, this%cp0, this%cum0, this%cm0
            close(io)
		end if

    end subroutine driven_current_result_save       
end module driven_current_module