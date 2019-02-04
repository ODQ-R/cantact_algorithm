PROGRAM contact_detection
	implicit none
	
	integer			:: i, j
	real			:: dij, distancia
	real			:: Eij
	real 			:: Rij
	real, parameter		:: E = 1e3
	real, parameter		:: poisson = 0.28
	integer, parameter	:: np = 4
	integer, parameter	:: c = 2
	integer, parameter	:: r = 3.00
	real, dimension (np,c)	:: xy
	real, dimension (np,c)	:: nij
	real, dimension (np,c)	:: Fij
	real, dimension (np,c)	:: Fji
	real, dimension (np,c)	:: V
	
	Eij = ((1 - poisson**2) / E) + ((1 - poisson**2) / E)
	Rij = (1.00 / r) + (1.00 / r)
	
	
	open (unit = 100, status = "old", file = 'C:\Users\osvaldo\Desktop\Códigos_2019\contact_detection\contact_detection\RANDOM_PARTICLE_XY.txt')
	open (unit = 200, status = "old", file = 'C:\Users\osvaldo\Desktop\Códigos_2019\contact_detection\contact_detection\RANDOM_PARTICLE_V.txt')
	
	do i = 1, np							
			read (100,*) (xy(i,j), j = 1, c)
			read (200,*) (V(i,j), j = 1, c)
			write (*,*)  (V(i,j), j = 1, c)
	end do
	
	
	do i = 1, np							
		do j = i+1, np				
			
			distancia = sqrt((xy(j,1) - xy(i,1))**2 + (xy(j,2) - xy(i,2))**2)
			dij = 2*r - distancia 
			
			nij(i,1) = ((xy(j,1) - (xy(i,1)))) / distancia
			nij(i,2) = ((xy(j,2) - (xy(i,2)))) / distancia
			
			if (dij > 0.00) then
				Fij(i,:) = (-4/3 * Eij * (Rij**2) * (dij)**3/2 ) * nij(i,:)
				Fij(i,:) = -Fij(i,:)
			else
				Fij(i,:) = 0.00
				Fij(i,:) = 0.00
			end if					
			
			V(i,:)  = V(i,:) !+ ((dt/m) * FT(i,1)) 
			
			!write(*,*) i,j,dij
		end do
	end do
	
	
	
	
	
	
	
close (100)
close (200)
END PROGRAM contact_detection

