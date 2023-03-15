module Utils
    use kind_module
    contains
    function sys_time()
    ! ** return system time
    implicit none
        real(wp) sys_time
        integer count, count_rate, count_max
        call system_clock(count, count_rate, count_max)
        sys_time = count*1.0/count_rate
    return
    end
end  module Utils