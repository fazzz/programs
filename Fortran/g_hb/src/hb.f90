module hydrogenbond
  implicit none

  contains

    integer function baker_hubbard ( pos_D, pos_H, pos_A, len, ang )
      implicit none
      
      real(8) :: pos_D(3), pos_H(3), pos_A(3)
      real(8) len, ang

      integer i, j, k

      len = calc_length(pos_H, pos_A)
      ang = calc_angle(pos_D, pos_A, pos_H)

      if (len < len_max .and. ang < ang_max) then
         baker_hubbard = 1
      else then
         baker_hubbard = 0
      end if
         
    end subroutine baker_hubbard
      
end module hydrogenbond
