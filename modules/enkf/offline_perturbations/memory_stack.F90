module memory_stack
!
! The memory_stack module provides a stack with numbers that can be used as index
! pointing to an array that stores some kind of data. Initially all elements are
! marked as free. Free elements can be requested for use (pop) and put back when
! no longer needed (push).
! The use of a stack with free memory slots avoids searching the storage array for
! empty space.
!

!
! The memstack type stores data for a memory stack.
! The main idea is to store to another array that is initially filled with
! the available indices. New values can be retrieved from the top or put back
! on the top
!

   implicit none

   type memstack
      integer, allocatable :: free_slots(:)
      integer :: current_size
      integer :: number_elements
   end type memstack

contains

   function memstack_create(initial_size)
      integer, intent(IN) :: initial_size
      type(memstack) :: memstack_create
      integer :: i
      memstack_create%current_size = initial_size
      memstack_create%number_elements = initial_size
      if (allocated(memstack_create%free_slots)) then
         deallocate (memstack_create%free_slots)
      end if
      allocate (memstack_create%free_slots(initial_size))
      do i = 1, memstack_create%current_size
         memstack_create%free_slots(i) = i
      end do

   end function

   function memstack_pop(stack)
      type(memstack), intent(INOUT) :: stack
      integer :: memstack_pop
      if (stack%current_size >= 1) then
         memstack_pop = stack%free_slots(stack%number_elements)
         stack%number_elements = stack%number_elements - 1; 
      else
         memstack_pop = -1
      end if
   end function

   subroutine memstack_push(stack, indx)
      type(memstack), intent(INOUT) :: stack
      integer, intent(IN) :: indx
      if (stack%number_elements >= stack%current_size) then
         write (*, *) 'Attempt to push more elements than fit on the stack'
         call exit(1)
      end if
      stack%number_elements = stack%number_elements + 1; 
      stack%free_slots(stack%number_elements) = indx
   end subroutine

   subroutine memstack_print(stack)
      type(memstack), intent(IN) :: stack
      integer :: i
      write (*, *) 'dump of memory_stack'
      write (*, *) '========================'
      write (*, *) 'size = ', stack%current_size
      write (*, *) 'number_elements = ', stack%number_elements
      write (*, *) '========================'
      do i = 1, stack%number_elements
         write (*, *) stack%free_slots(i)
      end do
      write (*, *) '========================'
   end subroutine

   subroutine memstack_free(stack)
      type(memstack), intent(INOUT) :: stack
      if (allocated(stack%free_slots)) then
         deallocate (stack%free_slots)
      end if
   end subroutine

end module
