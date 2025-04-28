!  reconstruction of water molecules          
         call reconstruction()

! Initialization of the reciprocal vectors   
         call init_fact_ew()  
! Calculation of the energy (vdw+coulombic) in the real space       (vdw+elec)   
	     call Inter_Realspace()
!  Calculation of the electrostatic energy in the reciprocal space   (Reci)  	     
         call Inter_Reci(xxx,yyy,zzz,xmolcm,ymolcm,zmolcm,0)
! Store the box length before the perturbation         
	     boxl_old(1)=boxl(1)  # length box x direction
	     boxl_old(5)=boxl(5)  # length box y direction
	     boxl_old(9)=boxl(9)  # length box z direction	     
!  vdw + electrostatic energy in the reference state 	     
	     U_LJ_0=Etot_tot_0  ! vdw+elec
	     Reci0=Reci_elec    ! Reci

! application of the Perturbation - Lx
             call trans_Area( X_new,Y_new,Z_new,Xmol_New,Ymol_New,Zmol_New,1)
! Same sequence that in the reference state             
             call init_fact_ew()
             call energy_TA(1,X_new,Y_new,Z_new,Xmol_New,Ymol_New,Zmol_New,U_LJ_1)
             call Inter_Reci(X_new,Y_new,Z_new,Xmol_New,Ymol_New,Zmol_New,1)
             boxl(1)=boxl_old(1)
             boxl(5)=boxl_old(5)
             boxl(9)=boxl_old(9)
             boxl9h=boxl(9)+diff
!  vdw + electrostatic energy in the perturbed state              
             Reci1=Reci_elec
             U_LJ_1=U_LJ_1             

!  epsilon_TA is the amplitude of the perturbation of Lx (box length according to the x direction)
             Delta_Lx=boxl(1)*epsilon_TA

!  average ensemble
             av_diffU2=av_diffU2 +exp(-(U_LJ_1 - U_LJ_0)/(8.314*0.001*temperature))
             dF_dx=-((8.314*temperature*0.001)/(Delta_Lx))*
     &log(((1+epsilon_TA)**(Nbmolec))*(av_diffU2/real(Nfile)))

!
!Nbmolec = number of molecules
!Nfile = number of configurations
!epsilon_TA=amplitude of the perturbation
!xxx,yyy,zzz - atomic coordinates in the reference state (denoted as 0)
!xmolcm,ymolcm,zmolcm - Coordinates of the center of mass in the reference state (labeled 0)
!boxl(1), boxl(5), boxl(9) - box lengths according to x, y and z directions respectively
!X_new,Y_new,Z_new  - atomic coordinates in the perturbed state (denoted as 1)
!Xmol_New,Ymol_New,Zmol_New - Coordinates of the center of mass in the perturbed state (labeled 1)

!Energy is in kJ/mol
!dF_dx in in kJ/mol/Angstrom

LAMMPS and MD-analysis tools are capable of providing van der Waals, real-space, and 
reciprocal-space electrostatic energies
 
subroutine  trans_Area(X_Ntmp,Y_Ntmp,Z_Ntmp,Xmol_Ntmp,Ymol_Ntmp,Zmol_Ntmp,delta)
 
      
      integer :: delta
      real(kind=8), dimension(maxatm)         :: X_Ntmp,Y_Ntmp,Z_Ntmp 
      real(kind=8), dimension(maxmol)         :: Xmol_Ntmp,Ymol_Ntmp,Zmol_Ntmp
      REAL(KIND=8) :: rappx,rappy,rappz,dx,dy,dz
      REAL(KIND=8) :: incre1,incre2,boxl_pl(10),increx,increy,increz
      
    
      
     
!Perturbation application       
         increx =  (1.0+epsilon_TA)
         increy =  1.
         increz =  1. 
        
!Effect on box length    
         boxl(1)=boxl_pl(1)* increx
         boxl(5)=boxl_pl(5)* increy
         boxl(9)=boxl_pl(9)* increz
        
!Effect on the center of mass position   
       compt=0
       do k=1,nmoltype
         do j=1,Nbmolec(k)
            compt=compt+1
            Xmol_Ntmp(compt)=xmolcm2(compt)*increx
	        Ymol_Ntmp(compt)=ymolcm2(compt)*increy
	        Zmol_Ntmp(compt)=zmolcm2(compt)*increz	       
	 enddo
       enddo       
 
!Atomistic positions resulting from the translation of the center of mass of molecules
      compt=0
      do k=1,nmoltype
         do j=1,Nbmolec(k)
            compt=compt+1
            do i =1,Lengt(k)
	           X_Ntmp(indice_atm(k,j,i))=xxx(indice_atm(k,j,i))+(Xmol_Ntmp(compt)-xmolcm2(compt))
               Y_Ntmp(indice_atm(k,j,i))=yyy(indice_atm(k,j,i))+(Ymol_Ntmp(compt)-ymolcm2(compt))
	           Z_Ntmp(indice_atm(k,j,i))=zzz(indice_atm(k,j,i))+(Zmol_Ntmp(compt)-zmolcm2(compt))
	      
	        enddo
	      enddo 
	  enddo	    
       
       
       
      
      end subroutine trans_Area
