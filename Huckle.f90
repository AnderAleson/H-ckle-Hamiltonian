program huckle

	!! VAIRABLE DECLARATION
	implicit none
	
	
	!!real*8 :: alpha
	integer :: num_atoms , num_e , charge, number_of_points,i, num_atom_type ,num_bond_type
	character(len=10) :: option,file_name
	real*8, dimension(:), allocatable :: eigen_values,beta,alpha
	!integer, dimension(:), allocatable :: title
	real*8, dimension(:,:), allocatable :: hamiltonian, wfn, wfnt

	!!END VARIABLE DECLARATION

	
	call get_command_argument(1, file_name)
	
	!!OPENING THE FILE IN READING MODE
	open(unit=1, file=file_name,  action="read")
	open(unit=2, file="output.txt",  action="write")
	open(unit=3, file="eigenvalues.out", action="write")
	open(unit=4, file="eigenvectors.out",  action="write")
	

	
	!!READ NUMBER OF ATOMS
	  READ(1,*)
	  READ(1,*)  num_atoms
	  write(*,*) "NUMBER OF CARBONS"
	  WRITE(*,*) num_atoms
	!!READ IF ITS AND OPEN OR CLOSED SYSTEM
	  READ(1,*)
	  READ(1,*) option
	  write(*,*)
	  WRITE(*,*) "THE SYSTEM IS AN: ",option
	!! READ THE NUMBER OF ATOM TYPE
	  READ(1,*)
	  READ(1,*) num_atom_type
	  write(*,*)
	  WRITE(*,*) "NUMBER OF ATOM TYPES"
	  write(*,*) num_atom_type
	  
	!!ALLOCATE ALPHA
	  allocate(alpha(num_atoms))
	  alpha(:)=0.d8
	  
	  do i=1, num_atom_type
	  	read(1,*)
	  	read(1,*) alpha(i)
	  end do
	  
	  i=num_atom_type+1
	  do while( i<= num_atoms)
	  	alpha(i)=alpha(i-num_atom_type)
	  	i=i+1
	  end do
	  
	  WRITE(*,*) "THE ATOM PARAMETERS ARE:"
	  WRITE(*,*) alpha(:)	
	!! READ THE number of bond type
	  READ(1,*)
	  READ(1,*) num_bond_type
	  write(*,*)
	  WRITE(*,*) "THE NUMBER OF BOND TYPES IS:"
	  write(*,*) num_bond_type
	  
	  !!ALLOCATE BETA
	  allocate(beta(num_atoms))
	  beta(:)=0.d8
	  
	  do i=1, num_bond_type
	  	read(1,*)
	  	read(1,*) beta(i)
	  end do
	  
	  i=num_bond_type+1
	  do while( i<= num_atoms)
	  	beta(i)=beta(i-num_bond_type)
	  	i=i+1
	  end do
	  
	  WRITE(*,*) "THE BOND PARAMETERS ARE:"
	  WRITE(*,*) beta(:)	  
	  
	  !!LETS START BUILDING THE HAMILTONIAN.
	  
	  !! SET THE DIMENSIONS OF THE HAMILTONIAN.
	  allocate(hamiltonian(num_atoms,num_atoms))
	  hamiltonian(:,:)=0.d8
	  
	  if (option=="closed") then
	  
	        do i=1, num_atoms
	  		hamiltonian(i,i)=alpha(i)
	  	end do
	        
	  	do i=1, num_atoms-1
	  		hamiltonian(i+1,i)=beta(i)
	  		hamiltonian(i,i+1)=beta(i)
	  	end do
	  	hamiltonian(1,num_atoms)=beta(i)
	  	hamiltonian(num_atoms,1)=beta(i)
	  	
	  	!!write hamiltonian to si if okey
	  	WRITE(*,*) "THE HUCKLE HAMILTONIAN OF THE SYSTEM"	
	  	do i=1, num_atoms
	  		write(*,101) hamiltonian(i,:)
	  	end do
	  	
	  
	  else if (option=="open") then
	  
	  	do i=1, num_atoms
	  		hamiltonian(i,i)=alpha(i)
	  	end do
	  	
	  	do i=1, num_atoms-1
	  		hamiltonian(i+1,i)=beta(i)
	  		hamiltonian(i,i+1)=beta(i)
	  	end do
	  	!!write hamiltonian to si if okey	
	  	WRITE(*,*) "THE HUCKLE HAMILTONIAN OF THE SYSTEM"
	  	do i=1, num_atoms
	  		write(*,101) hamiltonian(i,:)
	  	end do
	  
	  else 
	  	WRITE(*,*) "Error reading option"
	  	stop
	  end if
	  
	  
	  !!NOW THAT WE HAVE SET THE HAMILTONIAN
	  !!LETS DIAGONALIZE IT
	  
	  !!LETS SET THE DIMENSIONS FOR THE WFN
	  allocate(wfn(num_atoms,num_atoms))
	  allocate(wfnt(num_atoms,num_atoms))
	  allocate(eigen_values(num_atoms))
	  
	  
	  wfn(:,:)=hamiltonian(:,:)
	  
	  call diagonalize_matrix(num_atoms,wfn, eigen_values)
	  wfnt=transpose(wfn)
	  
	  WRITE(*,*) "WRITING THE EIGEN VALUES"
	  WRITE(2,*) "WRITING THE EIGEN VALUES"
	  do i=1 , num_atoms
		write(*,102) "eigenvalue ",i,eigen_values(i)
		write(2,102) "eigenvalue ",i,eigen_values(i)
		write(3,*) i, eigen_values(i)
	  end do 
	  
	  WRITE(*,*) "WRITING THE EIGEN VECTORS"
	  WRITE(2,*) "WRITING THE EIGEN VECTORS"
	  
	  
	  do i=1 , num_atoms
		write(*,104)"eigenvector ", i, wfnt(i,:)
		write(2,104)"eigenvector ", i, wfnt(i,:)
		!write(4,105) i,  wfn(i,:)
		write(4, '(I5, 5X, 1100F10.4)') i, wfn(i, :)
	  end do 
	  
	  call plot_MO_py(eigen_values,num_atoms,file_name)	  

	  
	  
	  

	close(1)
	close(2)
	close(3)
	close(4)
	stop
	!! define writing formats
101	format(4x,20f8.4)
105	format(i4,4x,20f8.4)
104     format(a13, i3,4x,1100f8.4)
102	format(a13,i10, 5x,f8.4)
!!103	format(i.3)



	end program huckle
	
	
	
	!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!!!!
subroutine diagonalize_matrix(N,A,e)

      ! Diagonalize a square matrix
    
      implicit none
    
      ! Input variables
    
      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N)
      double precision,intent(out)  :: e(N)
    
      ! Local variables
    
      integer                       :: lwork,info
      integer                       :: i
      double precision,allocatable  :: work(:)
      
    
      ! Memory allocation
    
      allocate(work(3*N))
      lwork = size(work)
    
      call dsyev('V','U',N,A,N,e,work,lwork,info)
     
      if(info /= 0) then 
        write(*,'(a)') 'Problem in diagonalize_matrix (dsyev)!!'
        stop
      endif
      
      do i = 1 , N
        if (abs(e(i)) < 1e-10) e(i) = 0
      end do  


	
	
	
end subroutine diagonalize_matrix

subroutine plot_MO_py(e,N,file_name)
   implicit none
   integer,intent(in) :: N
   double precision,intent(in):: e(N)
   character(len=10),intent(in) :: file_name
   
   integer :: i,length
   character(len=10) :: output_filename
   
   length = len_trim(file_name)
   
   do i = length, 1, -1
        if (file_name(i:i) == '.') then
            output_filename = file_name(1:i-1) // '.py'
            exit
        end if
    end do
   

   open(unit=19, file=output_filename,  action="write")
   
   write(19,'(A)') "import matplotlib.pyplot as plt"
   write(19,'(A)') "import numpy as np"
   write(19,'(A)') "from collections import Counter"
   write(19,'(A)') "energy_levels = np.array(["
   do i=1,N-1
      write(19,*) e(i),","
   end do
   write(19,*) e(N),"])"
   write(19,'(A)') "x0 = 0  "
write(19,'(A)') "dx = 0.3  "

write(19,'(A)') "x_positions = []"
write(19,'(A)') "unique_energy_levels = []"

write(19,'(A)') "i = 0"
write(19,'(A)') "while i < len(energy_levels):"
write(19,'(A)') "    if i < len(energy_levels) - 1 and abs(energy_levels[i] - energy_levels[i + 1]) < 10**-8:"
write(19,'(A)') "        # Degenerate case: split positions"
write(19,'(A)') "        x_positions.append(x0 - dx)"
write(19,'(A)') "        x_positions.append(x0 + dx)"
write(19,'(A)') "        unique_energy_levels.append(energy_levels[i])"
write(19,'(A)') "        unique_energy_levels.append(energy_levels[i + 1])"
write(19,'(A)') "        i += 2  # Skip next value since it's already plotted"
write(19,'(A)') "    else:"
write(19,'(A)') "        # Non-degenerate case: plot at center"
write(19,'(A)') "        x_positions.append(x0)"
write(19,'(A)') "        unique_energy_levels.append(energy_levels[i])"
write(19,'(A)') "        i += 1"
write(19,'(A)') "    colors = ['b' if e < 0 else 'g' for e in unique_energy_levels]"

write(19,'(A)') ""
write(19,'(A)') "# Plotting"
write(19,'(A)') "plt.figure(figsize=(5, 6))"
write(19,'(A)') "plt.scatter(x_positions, unique_energy_levels, color=colors, marker='_', s=200)"
write(19,'(A)') "plt.xlabel('')"
write(19,'(A)') "plt.ylabel('Energy [Hartree]')"
write(19,'(A)') "plt.title('MO diagram')"
write(19,'(A)') "plt.grid(True, linestyle='--', alpha=0.6)"
write(19,'(A)') "plt.show()"

 close(19)


end subroutine plot_MO_py
	
	!!!!!!!!!! END SUBROUTINES !!!!!!!!!!!!


