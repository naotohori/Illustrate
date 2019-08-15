!                 ******************************
!                      I L L U S T R A T E
!                   Biomolecular Illustration
!       ******************************
!       copyright 2019 David S Goodsell
!
! Licensed under the Apache License, Version 2.0 (the "License");  
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License
!
! This work was supported by Damon Runyon-Walter Winchell Cancer Research Fund Fellowship DRG 972,
! the US National Institutes of Health R01-GM120604 and the kind support of the RCSB Protein Data Bank.
!
!--------------------------------------------------------------------
!
! DS Goodsell & AJ Olson (1992) "Molecular Illustration in Black and White" JMolGraphics 10, 235-240.
! April 2019 -- Simplified and released with only non-photorealistic rendering
!
!--------------------------------------------------------------------
! compile: gfortran illustrate.f -o illustrate
! run: illustrate < command_file
!
! Command file is read from unit 5
! Command file has command cards, followed by parameter cards
! Idiosyncracies (warning, postdoc code):
! *** must issue command cards in this order
! *** any number of rotation cards may be added, and they are concatenated when added
!      this means they are effectively applied last to first
!      so if you're progressively refining a position, add new rotations to the top of the list
! *** origin at upper left, +x down, +y left to right, +z towards viewer, molecules clipped at z=0
!
! read                                          #READ command
! 2hhb.pdb                                      #PDB format coordinate file
! HETATM-----HOH-- 0,9999, 0.5,0.5,0.5, 0.0     #selection/rendering cards
! ATOM  -H-------- 0,9999, 0.5,0.5,0.5, 0.0
! ATOM  H--------- 0,9999, 0.5,0.5,0.5, 0.0
! ATOM  -C-------A 0,9999, 1.0,0.6,0.6, 1.6
! ATOM  -S-------A 0,9999, 1.0,0.5,0.5, 1.8
! ATOM  ---------A 0,9999, 1.0,0.5,0.5, 1.5
! ATOM  -C-------C 0,9999, 1.0,0.6,0.6, 1.6
! ATOM  -S-------C 0,9999, 1.0,0.5,0.5, 1.8
! ATOM  ---------C 0,9999, 1.0,0.5,0.5, 1.5
! ATOM  -C-------- 0,9999, 1.0,0.8,0.6, 1.6
! ATOM  -S-------- 0,9999, 1.0,0.7,0.5, 1.8
! ATOM  ---------- 0,9999, 1.0,0.7,0.5, 1.5
! HETATMFE---HEM-- 0,9999, 1.0,0.8,0.0, 1.8
! HETATM-C---HEM-- 0,9999, 1.0,0.3,0.3, 1.6
! HETATM-----HEM-- 0,9999, 1.0,0.1,0.1, 1.5
! END                                          #end of READ commad
! center                                       #CENTER command 
! auto                                         #use autocentering (typical)
! trans                                        #TRANSLATION command
! 0.,0.,0.                                     #x,y,z for translation
! scale                                        #SCALE command
! 12.0                                         #scale value (pixels/Angstrom)
! zrot                                         #ROTATION command
! 90.                                          #rotation angle (deg)
! wor                                          #WORLD rendering parameter
! 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0              #rgb background, rgb fog, fractional opacity of fog front & back
! 1,0.0023,2.0,1.0,0.2                         #soft shadow parameters
! -30,-30                                      #image size, negative for autosizing with padding
! illustrate                                   #ILLUSTRATION command
! 3.0,10.0,4,0.0,5.0                           #contour outlines
! 3.0,10.0                                     #subunit outlines
! 3.0,8.0,6000.                                #residue outlines
! calculate                                    #CALCULATE command
! 2hhb.pnm                                     #image file name in ppm format
!
!--------------------------------------------------------------------
   integer*4 iost
   logical flg_cycle
! ***** OUTPUT SCAN LINE *****
   integer*4 scanline(9000)
! ***** FRAME BUFFER WITH COLORS *****
   real*4 pix(-10:3008,-10:3008,4)
! ***** FRAME BUFFER CONTAINING THE z HEIGHT *****
   real*4 zpix(-10:3008,-10:3008)
! ***** FRAME BUFFER CONTAINING THE ATOM NUMBER *****
   integer*4 atom(-10:3008,-10:3008)
   integer*4 bio(-10:3008,-10:3008)
   integer*4 n,ia
! ***** SHADING MAPS FOR THE EIGHT ATOM TYPES *****
   real*4 sphdat(0:32000,3)
   integer*4 numpix
! ***** PHONG SHADING PARAMETERS ****** 
   real*4 colortype(0:1000,3)
   real*4 rback(3),rfog(3)
! ***** ATOMIC INFORMATION *****
   real*4 coord(3,350000)
   integer*4 type(350000),res(0:350000),su(0:350000)
   real*4 radtype(1000)
   character(1) chain,chainlast
! ***** transformation matrices *****
   real*4 matrixin(4,4),rm(4,4)
! ***** biological unit stuff *****
   real*4 biomat(4,4,500)
   integer*4 nbiochain,nbiomat
   character(1) biochain(500)
! ***** STUFF FOR OUTLINES *****
   real*4 l_opacity,l_opacity_ave,g_opacity,opacity
   real*4 l(-1:1,-1:1),g
   real*4 l_low,l_high,g_low,g_high
   real*4 l_diff_min, l_diff_max
! ***** Conical shadows *****
   real*4 rtable(-51:51,-51:51),coneangle,pcone,rcone
! ***** ETC. *****
   real*4 x,y,z,rx,ry,rz,xp,yp,zp,d
   real*4 xn,yn,zn,ci,xs,xy,xz,cs
   character(3) comcode(10),command
   character(80) filename,inputfile
   integer*4 ixsize,iysize, idepth
   integer*4 ix,iy,iz
   integer illustrationflag
! ***** stuff for reading atoms
   integer*4 resrange(2,1000),inptype(1000),inpsu(1000)
   character(6) atomdescriptor(1000)
   character(10) descriptor(1000)
   character(80) instring
!--------------------------------------------------------------------
! ***** Initialize a few things *****
!--------------------------------------------------------------------
   rscale=1.

   illustrationflag=0

   comcode = ['rea','tra','xro','yro','zro','sca','cen','wor','cal','ill']

   call clearmatrix(rm)
   l_low=1.
   l_high=10.
   g_low=350000.
   g_high=21000.
   xshmax=4000
   yshmax=4000
   su(0)=9999
   res(0)=9999
   xtran=0.
   ytran=0.
   ztran=0.
   l_diff_max=50.
   l_diff_min=1.
   ixsize=0
   iysize=0
   do i=1,3
      do j=1,4
         do k=1,500
            biomat(i,j,k)=0.
            if (i.eq.j) biomat(i,j,k)=1.
         enddo
      enddo
   enddo

! ********************* READ CONTROL CARDS ***************************
   do
      read (5,'(a3)',iostat=iost) command
      if (iost < 0) exit

      icommand=10
      do icount=1,10
         if (command.eq.comcode(icount)) icommand=icount
      enddo

      select case (icommand)
   
      case default
         write (6,'(a35,a3,a7)') ' ***** invalid control card read: ', command, ' ***** '

      !--------------------------------------------------------------------
      !  *** read and classify atoms ***
      !
      case (1)
         read(5,'(a80)') inputfile
         open(1,file=inputfile,form='formatted',status='old')
         ! --- read atom descriptors ---
         ! param: atom descriptor cards
         ! ATOM  -CA--ALA-A 0,9999,0.5,0.5,0.5,1.6
         ! ^^^^^^   ATOM or HETATM, matched with columns 1-6 of PDB record
         !       ^^^^^^^^^^ compared to columns 13-22 of PDB record
         !                  "-" is a wildcard
         !                  ^^^^^^ residue number range
         !                         ^^^^^^^^^^^ r,g,b (0.-1.)
         !                                     ^^^ radius (Angstrom)
         ! cards are read in order, and if there is a match, the atom is assigned that type
         ! cards are read until one is found without "ATOM  " or "HETATM"
         !
         ndes=0
         ! --- default 50% gray ---
         do i=0,1000
            do j=1,3
               colortype(i,j)=.5
            enddo
         enddo
         do
            read(5,'(a80)') instring
            if (instring(1:3).eq.'END') exit
            ndes=ndes+1
            read(instring,'(a6,a10)') atomdescriptor(ndes),descriptor(ndes)
            read(instring(18:80),*) (resrange(i,ndes),i=1,2),rr,rg,rb,rad
            colortype(ndes,1)=rr
            colortype(ndes,2)=rg
            colortype(ndes,3)=rb
            radtype(ndes)=rad
         enddo
         write(6,*) ' atom descriptors: ',ndes
         do i=1,ndes
            write(6,*) "type, color, radius ",i,(colortype(i,j),j=1,3),radtype(i)
         enddo
   
         ! --- read atoms and classify ---
         n=0
         nsu=0
         chain=" "
         nbiomat=0
         nbiochain=0
   
         do
            read(1,'(a80)',iostat=iost) instring
            if (iost < 0) exit
   
            if (instring(1:5).eq."MODEL") nsu=nsu+1
            if (instring(12:25).eq."BIOMOLECULE: 1") then
               do
                  read(1,'(a80)') instring
                  if (instring(12:16).eq."APPLY") then
                     ich=43
                     do while (instring(ich:ich).ne." ")
                        nbiochain=nbiochain+1
                        biochain(nbiochain)=instring(ich:ich)
                        ich=ich+3
                     enddo
                     write(6,*) "Chains in Biological Assembly"
                     do i=1,nbiochain
                        write(6,*) i,biochain(i)
                     enddo
                     write(6,*)
                  endif
                  if (instring(14:19).eq."BIOMT1") then
                     nbiomat=nbiomat+1
                     read(instring(20:80),*) ib,(biomat(1,j,nbiomat),j=1,4)
                  endif
                  if (instring(14:19).eq."BIOMT2") then
                     !read(instring,'(23x,5f10.6,f15.5)') (biomat(2,j,nbiomat),j=1,4)
                     read(instring(20:80),*) ib,(biomat(2,j,nbiomat),j=1,4)
                  endif
                  if (instring(14:19).eq."BIOMT3") then
                     !read(instring,'(23x,5f10.6,f15.5)') (biomat(3,j,nbiomat),j=1,4)
                     read(instring(20:80),*) ib,(biomat(3,j,nbiomat),j=1,4)
                  endif
                  if (instring(14:19).eq."      ") then
                     write(6,*) "Number of BIOMT ",nbiomat
                     do ibio=1,nbiomat
                        do im=1,3
                           write(6,*) "BIOMAT ",ibio,(biomat(im,j,ibio),j=1,4)
                        enddo
                     enddo
            
                     exit
                  endif
               enddo ! --- end of BIOMAT
            endif
   
            if ((instring(1:4).ne.'ATOM').and.(instring(1:6).ne.'HETATM')) cycle
   
            read(instring,'(22x,i4)') ires
   
            do ides=1,ndes
   
               if (instring(1:6).ne.atomdescriptor(ides)(1:6)) cycle
   
               flg_cycle = .false.
               do ia=1,10
                  if (descriptor(ides)(ia:ia).eq.'-') cycle
                  if (instring(12+ia:12+ia).ne.descriptor(ides)(ia:ia)) then
                     flg_cycle = .true.
                     exit
                  endif
               enddo
               if (flg_cycle) cycle
      
               if ((ires.lt.resrange(1,ides)).or.(ires.gt.resrange(2,ides))) cycle
      
               if (radtype(ides).eq.0) exit
      
               ! --- check if chain is in biological assembly
               if (nbiochain.ne.0) then
                  ibioflag=0
                  do ibio=1,nbiochain
                     if (instring(22:22).eq.biochain(ibio)) ibioflag=1
                  enddo
                  if (ibioflag.eq.0) exit
               endif
      
               ! --- found an atom to save ---
               n=n+1
               read(instring,'(30x,3f8.3)') (coord(i,n),i=1,3)
               type(n)=ides
               ! --- assign subunits automatically ---
               chain=instring(22:22)
               if (chain.ne.chainlast) then
                  nsu=nsu+1
                  chainlast=chain
               endif
               su(n)=nsu
               res(n)=ires
   
               ! --- apply biomat ---
               !if (nbiomat.gt.0) then
               !rx=coord(1,n)
               !ry=coord(2,n)
               !rz=coord(3,n)
               !do imat=1,nbiomat
               !su(n+imat)=su(n)
               !res(n+imat)=res(n)
               !type(n+imat)=ides
               !coord(1,n+imat)=rx*biomat(1,1,imat)+
               !    &                  ry*biomat(1,2,imat)+
               !    &                  rz*biomat(1,3,imat)+
               !    &                     biomat(1,4,imat)
               !coord(2,n+imat)=rx*biomat(2,1,imat)+
               !    &                  ry*biomat(2,2,imat)+
               !    &                  rz*biomat(2,3,imat)+
               !    &                     biomat(2,4,imat)
               !coord(3,n+imat)=rx*biomat(3,1,imat)+
               !    &                  ry*biomat(3,2,imat)+
               !    &                  rz*biomat(3,3,imat)+
               !    &                     biomat(3,4,imat)
               !enddo
               !n=n+nbiomat
               !endif
      
               exit
            enddo
         enddo 
         ! --- done reading atoms ---         
   
         write(6,*)' atoms read: ', n, ' from: ',inputfile
         write(6,*) " number of subunits: ",nsu
         write(6,*)' '
   
   !--------------------------------------------------------------------
   !       TRANSLATION 
      case (2)
         read(5,*) xtrani,ytrani,ztrani
         xtran=xtran+xtrani
         ytran=ytran+ytrani
         ztran=ztran+ztrani
         write(6,*) 'translation vector : ',xtran,ytran,ztran
         write(6,*)
   
      !--------------------------------------------------------------------
      !  Z ROTATION
      case (5)
         call clearmatrix(matrixin)
         read(5,*) angle
         write(6,*) 'z rotation : ',angle
         ! --minus sign is because original rotations were left-handed! (warning: postdoc code)
         angle=-angle*3.141592/180.
         matrixin(1,1)=cos(angle)
         matrixin(1,2)=-sin(angle)
         matrixin(2,1)=sin(angle)
         matrixin(2,2)=cos(angle)
         call catenate(rm,matrixin)
   
      !--------------------------------------------------------------------
      !  Y ROTATION
      case (4)
         call clearmatrix(matrixin)
         read(5,*) angle
         write(6,*) 'y rotation : ',angle
         angle=-angle*3.141592/180.
         matrixin(1,1)=cos(angle)
         matrixin(1,3)=sin(angle)
         matrixin(3,1)=-sin(angle)
         matrixin(3,3)=cos(angle)
         call catenate(rm,matrixin)
   
      !--------------------------------------------------------------------
      !  X ROTATION
      case (3)
         call clearmatrix(matrixin)
         read(5,*) angle
         write(6,*) 'x rotation : ',angle
         angle=-angle*3.141592/180.
         matrixin(2,2)=cos(angle)
         matrixin(2,3)=-sin(angle)
         matrixin(3,2)=sin(angle)
         matrixin(3,3)=cos(angle)
         call catenate(rm,matrixin)
   
      !--------------------------------------------------------------------
      !       SCALE 
      case (6)
         read(5,*) rscalei
         rscale=rscale*rscalei
         write(6,*) 'scale factor : ',rscale
         write(6,*)
   
      !--------------------------------------------------------------------
      !  CENTERING of coordinates
      !       three options: CEN(TER) will center on max/min coordinates in x,y,z
      !                      AUT(O-CENTER) will center the rotated coordinates
      !                        in the frame, and place the uppermost atom
      !                        just below the image plane
      !              (done when image processing begins)
      !       ROTATIONs are applied before the centering
      !       TRANSLATIONs are applied after the centering
      case (7)
         read(5,'(a3)') command
         autocenter=0
         if (command.eq.'aut') then
            autocenter=1
         else if (command.eq.'cen') then
            autocenter=2
         endif
      
      !--------------------------------------------------------------------
      !  WORLD parameters that describe the environment 
      case (8)
         ! param: rback(3) -- rgb of background, 0.-1.
         ! param: rfog(3) -- rgb of fog, 0.-1.
         ! param: pfogh,pfogl -- fractional intensity of fog color at front and back of molecule
         read (5,*) (rback(i),i=1,3),(rfog(i),i=1,3),pfogh,pfogl
         if (pfogl.lt.0.) then
            pfogl=-pfogl
         endif
      
         do i=1,3
         if (rback(i).gt.1.) rback(i)=1.
         if (rfog(i).gt.1.) rfog(i)=1.
              colortype(0,i)=rback(i)
         enddo
      
         if (pfogh.gt.1.) pfogh=1.
         if (pfogl.gt.1.) pfogl=1.
      
         write (6,'(1x,a20,1x,8f7.2)') ' background inten. :',(rback(i),i=1,3)
         write (6,'(1x,a20,1x,8f7.2)') ' fog intensity :    ',(rfog(i),i=1,3)
         write (6,'(1x,a20,1x,8f7.2)') ' upper fog percent :',pfogh*100.
         write (6,'(1x,a20,1x,8f7.2)') ' lower fog percent :',pfogl*100.
         pfogdiff=pfogh-pfogl
      
         ! param: icone -- flag, 0=no soft shadow, other=soft shadow
         ! param: pcone -- fractional shadowing contribution from each atom (0.0023)
         !                 larger=darker shadow
         ! param: coneangle -- angle of shadowing around each atom (2.0)
         !                     larger values give tighter shadowed region
         ! param: rcone -- shadowing applied only if z values of atoms are
         !                 greater than rcone (removes many small shadows in
         !                 creases) (1.0 A)
         ! param: pshadowmax -- maximal shadowing amount (0.7) smaller=darker shadow
         read (5,*) icone,pcone,coneangle,rcone,pshadowmax
         if (icone.ne.0) write (6,'(1x,a20,1x,8f7.2)') ' draw conical shadows'
      
         ! param: ixsize,iysize -- vertical and horizonal size (pixels)
         !        negative values are autosized, padded by the value 
         read(5,*) ixsize,iysize
         if (ixsize.gt.3000) ixsize=3000
         if (iysize.gt.3000) iysize=3000
         write(6,*) 'input value for image size', ixsize,iysize
   
      !--------------------------------------------------------------------
      !  ILLUSTRATION parameters for outlines
      case (10)
         illustrationflag=1
      
         ! param: l_low, l_high -- thresholds for drawing contour outlines (3.,8.)
         !                    higher values give fewer outlines
         !                         narrower range gives jaggier outlines
         !                    values are based on the number of pixels in the kernel
         ! param: ikernel -- kernel for doing derivative, values=1,2,3,4 (4)
         ! param: l_diff_min, l_diff_max -- range of difference in z-values used for calculation (0.,5. A)
         !                    0-1 gives outlines around every atom, 0-1000 outlines the whole molecule
         ! param: r_low, r_high -- thresholds for drawing subunit outlines (3.,8.)
         ! param: g_low, g_high -- thresholds for drawing residue outlines (3.,8.)
         ! param: resdiff -- residue outlines drawn for atoms with this difference in residue numbers or greater
      
         read(5,*) l_low,l_high,ikernel,l_diff_min,l_diff_max
         read(5,*) r_low,r_high
         read(5,*) g_low,g_high,resdiff
         write(6,*) 'illustration parameters'
         write(6,*) 'l parameters: ',l_low,l_high
         write(6,*) 'g parameters: ',g_low,g_high
   
      !--------------------------------------------------------------------
      case (9)
         write (6,'(1x,a20)') ' *begin calculation*'
   
         ! --- if no BIOMT in file, use biomat 1 == identity matrix
         nbiomat=max(nbiomat,1)
         ! ***** Populate conical shadow table ****
         conemax=50.
         do i=-51,51
            do j=-51,51
               rtable(i,j)=sqrt(float(i)**2+float(j)**2)
               if (rtable(i,j).gt.conemax) rtable(i,j)=10000.
            enddo
         enddo
         rtable(0,0)=10000.
         ! ***** SCALE RADII *****
         do i=1,ndes
            radtype(i)=radtype(i)*rscale
            if (radtype(i).gt.radius_max) radius_max=radtype(i)
         enddo
         ! ***** APPLY AUTOCENTERING and AUTOSIZING, if switched on *****
         if (autocenter.gt.0) then
   
            xmin=10000.
            xmax=-10000.
            ymin=10000.
            ymax=-10000.
            zmin=10000.
            zmax=-10000.
      
            do ia=1,n
               do ibio=1,nbiomat
                  !--apply biomat
                  rx=coord(1,ia)*biomat(1,1,ibio)+coord(2,ia)*biomat(1,2,ibio)+ &
                     coord(3,ia)*biomat(1,3,ibio)+biomat(1,4,ibio)
                  ry=coord(1,ia)*biomat(2,1,ibio)+coord(2,ia)*biomat(2,2,ibio)+ &
                     coord(3,ia)*biomat(2,3,ibio)+biomat(2,4,ibio)
                  rz=coord(1,ia)*biomat(3,1,ibio)+coord(2,ia)*biomat(3,2,ibio)+ &
                     coord(3,ia)*biomat(3,3,ibio)+biomat(3,4,ibio)
            
                  !--apply rotation matrix
                  rx2=rx*rm(1,1)+ry*rm(2,1)+rz*rm(3,1)
                  ry2=rx*rm(1,2)+ry*rm(2,2)+rz*rm(3,2)
                  rz2=rx*rm(1,3)+ry*rm(2,3)+rz*rm(3,3)
                  xmin=MIN(xmin,rx2)
                  xmax=MAX(xmax,rx2)
                  ymin=MIN(ymin,ry2)
                  ymax=MAX(ymax,ry2)
                  zmin=MIN(zmin,rz2)
                  zmax=MAX(zmax,rz2)
      
               enddo
            enddo
      
            write(6,*) 'min coordinates : ',xmin,ymin,zmin
            write(6,*) 'max coordinates : ',xmax,ymax,zmax
            xtranc=-xmin-(xmax-xmin)/2.
            ytranc=-ymin-(ymax-ymin)/2.
            if (autocenter.eq.1) then
               ztranc=-zmax-radius_max-1.
               write(6,*) 'automating centering'
            endif
            if (autocenter.eq.2) then
               ztranc=-zmin-(zmax-zmin)/2.
               write(6,*) 'x,y,z centering'
            endif
         
            write(6,*) 'centering vector : ',xtranc,ytranc,ztranc
         
            if ((ixsize.le.0).or.(iysize.le.0)) then
              write(6,*)
              write(6,*) 'applying autosizing'
              write(6,*) 'x and y frame width: ',-ixsize,-iysize
              ixsize=-2.*ixsize+2.*radius_max+(xmax-xmin)*rscale
              iysize=-2.*iysize+2.*radius_max+(ymax-ymin)*rscale
            endif
         endif
   
         ixsize=min(ixsize,3000)
         iysize=min(iysize,3000)
   
         ! ***** OPEN OUTPUT FILES *****
         ixsize=int(ixsize/2)*2
         iysize=int(iysize/2)*2
         write(6,*) 'xsize and ysize: ',ixsize,iysize
         write(6,*)
         read(5,'(a80)') filename
         write(6,*) "output pnm filename: ",filename
         open(8,file=filename,form='formatted')
         write(8,'(a2)') "P3"
         write(8,'(2i5)') iysize,ixsize
         write(8,'(2i5)') 255
         ! this will write an file with opacities
         ! open(9,file="opacity.pnm",form='formatted')
         ! write(9,'(a2)') "P3"
         ! write(9,'(2i5)') iysize,ixsize
         ! write(9,'(2i5)') 255
         ! ***** MAP SPHERICAL SURFACES OVER ATOMS *****
         do ix=1,ixsize
            do iy=1,iysize
               pix(ix,iy,1)=0.
               pix(ix,iy,2)=0.
               pix(ix,iy,3)=0.
               pix(ix,iy,4)=0.
               atom(ix,iy)=0
               bio(ix,iy)=1
               zpix(ix,iy)=-10000.
            enddo
         enddo
   
         if (n.gt.0) then
   
            ! ----- create the spherical shading map for atom types -----
            do irad=1,ndes
               ic=1
               irlim=int(radtype(irad))
               if (irlim > 100) then
                  write(6,*) 'atoms radius * scale > 100'
                  stop
               endif
   
               do ix=-irlim-1,irlim+1
                  do iy=-irlim-1,irlim+1
                     x=float(ix)
                     y=float(iy)
                     d=sqrt(x*x+y*y)
                     if (d > radtype(irad)) cycle
                     z=sqrt(radtype(irad)**2-d*d)
                     sphdat(ic,1)=x
                     sphdat(ic,2)=y
                     sphdat(ic,3)=z
                     ic=ic+1
                  enddo
               enddo
               numpix=ic-1
               ! ----- then map spherical surface over atoms of the proper type ------
               icount=0
   
               do ia=1,n
   
                  if (type(ia).ne.irad) cycle
   
                  icount=icount+1
                  do ibio=1,nbiomat
                  !--apply biomat
                     rx=coord(1,ia)*biomat(1,1,ibio)+coord(2,ia)*biomat(1,2,ibio)+&
                             coord(3,ia)*biomat(1,3,ibio)+biomat(1,4,ibio)
                     ry=coord(1,ia)*biomat(2,1,ibio)+coord(2,ia)*biomat(2,2,ibio)+&
                             coord(3,ia)*biomat(2,3,ibio)+biomat(2,4,ibio)
                     rz=coord(1,ia)*biomat(3,1,ibio)+coord(2,ia)*biomat(3,2,ibio)+&
                             coord(3,ia)*biomat(3,3,ibio)+biomat(3,4,ibio)
      
                     !--apply rotation matrix
                     rx2=rx*rm(1,1)+ry*rm(2,1)+rz*rm(3,1)
                     ry2=rx*rm(1,2)+ry*rm(2,2)+rz*rm(3,2)
                     rz2=rx*rm(1,3)+ry*rm(2,3)+rz*rm(3,3)
                  
                     !--apply centering vector
                     rx2=rx2+xtranc
                     ry2=ry2+ytranc
                     rz2=rz2+ztranc
                  
                     !--apply translation and scaling
                     rx2=(rx2+xtran)*rscale
                     ry2=(ry2+ytran)*rscale
                     rz2=(rz2+ztran)*rscale
   
                     if (rz2 < 0.) then
                        do ipix=1,numpix
                           x=sphdat(ipix,1)+rx2+float(ixsize)/2.     
                           y=sphdat(ipix,2)+ry2+float(iysize)/2.
                           ix=int(x)
                           iy=int(y)
                           if ((x.gt.float(ixsize)).or.(x.lt.1.).or.&
                               (y.gt.float(iysize)).or.(y.lt.1.)) cycle
                           z=sphdat(ipix,3)+rz2
                           if (z.gt.zpix(ix,iy)) then
                              zpix(ix,iy)=z
                              atom(ix,iy)=ia
                              bio(ix,iy)=ibio
                           endif
                        enddo
                     endif
   
                  enddo ! -- ibio loop
               enddo ! -- ia loop
   
               write(6,*) icount,' spheres added of type: ',irad
            enddo ! -- irad loop
   
            write(6,*) 'shading maps written into depth buffer'
         endif
   
         !***** CALCULATE SECOND DERIVATIVE OUTLINES ******
         ! ---- find maximum and minimum z levels ---
         !  (note: I use a value of zpix=-10000. to distinguish background)
         zpix_max=-100000.
         zpix_min=100000.
         do ix=1,ixsize
            do iy=1,iysize
               if (zpix(ix,iy).gt.zpix_max) zpix_max=zpix(ix,iy)
               if ((zpix(ix,iy).ne.-10000.).and.&
                   (zpix(ix,iy).lt.zpix_min)) zpix_min=zpix(ix,iy)
            enddo
         enddo
         zpix_max=min(zpix_max,0.)
         zpix_spread=zpix_max-zpix_min
         write(6,*) 'zpix_min,zpix_max ',zpix_min,zpix_max
   
         ! ***** PROCESSING OF THE IMAGE BEGINS HERE*****
         l_diff_min=l_diff_min*rscale
         l_diff_max=l_diff_max*rscale
         write(6,*) ' Pixel processing beginning '
         do ix=1,ixsize
            do iy=1,iysize
               zpix(ix,iy)=min(zpix(ix,iy),0.)
            enddo
         enddo
   
         do ix=1,ixsize
            do iy=1,iysize
               ! ***** CONICAL SHADOW TESTING *****
               pconetot=1.
               if ((icone.ne.0).and.(atom(ix,iy).ne.0)) then   
                  do i=-50,50,5
                     do j=-50,50,5
                        if ((ix+i.gt.0).and.(ix+i.lt.ixsize).and.&
                            (iy+j.gt.0).and.(iy+j.lt.iysize)) then
                           rzdiff=zpix(ix+i,iy+j)-zpix(ix,iy)
                           if (rzdiff > rcone) then
                              if (rtable(i,j)*coneangle < rzdiff+rcone) then
                                 pconetot=pconetot-pcone
                              endif
                           endif
                        endif
                     enddo
                  enddo
                  pconetot=max(pconetot,pshadowmax)
               endif
               ! ***** CALCULATE THE FOG PERCENTAGE (DEPTH CUEING) *****
               pfh=pfogh-(zpix_max-zpix(ix,iy))/zpix_spread*pfogdiff
               if (zpix(ix,iy) < zpix_min) pfh=1.
               ! ***** CALCULATE OUTLINES *****
               g_opacity=0.
               l_opacity=0.
               if (illustrationflag /= 0) then
                  ! ***** SUBUNIT OUTLINES *****
                  ! ---- (this replaces original calculation of first derivatives)
   
                  if ((ix > 1).and.(ix < ixsize).and.(iy > 1).and.(iy < iysize)) then
                     ! --- calculate subunit 'derivatives'  ---
                     g=0.
                     r=0.
                     do i=-2,2
                        do j=-2,2
                           if (abs(i*j).ne.4) then
                              if ((su(atom(ix,iy)).ne.su(atom(ix+i,iy+j))).or.&
                                  (bio(ix,iy).ne.bio(ix+i,iy+j))) r=r+1.
                              if (abs(res(atom(ix,iy))-res(atom(ix+i,iy+j))) > resdiff) g=g+1.
                           endif
                        enddo
                     enddo
                     ! --- opacities are 1 for completely opaque ---
                     g_opacity=min((g-g_low)/(g_high-g_low),1.)
                     r_opacity=min((r-r_low)/(r_high-r_low),1.)
                     g_opacity=max(g_opacity,r_opacity,0.)
                  endif
                  ! ***** SECOND DERIVATIVE OUTLINES *****
                  if ((ix.gt.2).and.(ix.lt.ixsize-1).and.&
                      (iy.gt.2).and.(iy.lt.iysize-1)) then
                     rl=0.
                     l_opacity_ave=0.
                     do ixl=-1,1
                        do iyl=-1,1
                           ixc=ix+ixl
                           iyc=iy+iyl
   
                           if (ikernel.eq.1) then
                              l(ixl,iyl)=abs(1./3. * ( &
                                  -0.8*zpix(ixc-1,iyc-1)-1.*zpix(ixc-1,iyc)-0.8*zpix(ixc-1,iyc+1)-&
                                  1.0* zpix(ixc,iyc-1)+7.2*zpix(ixc,iyc)-1.0*zpix(ixc,iyc+1)-     &
                                  0.8* zpix(ixc+1,iyc-1)-1.*zpix(ixc+1,iyc)-0.8*zpix(ixc+1,iyc+1) &
                                   ))                                                             
                           endif
   
                           if (ikernel.eq.2) then
                              l(ixl,iyl)=abs(1./3. * ( &
                                  -0.8*zpix(ixc-1,iyc-1)-1.0*zpix(ixc-1,iyc)-0.8*zpix(ixc-1,iyc+1)-&
                                   1.0*zpix(ixc,iyc-1)+8.8*zpix(ixc,iyc)-1.0*zpix(ixc,iyc+1)-&
                                              0.8*zpix(ixc+1,iyc-1)-1.0*zpix(ixc+1,iyc)-0.8*zpix(ixc+1,iyc+1)-&
                                   0.1*zpix(ixc+2,iyc-1)-0.2*zpix(ixc+2,iyc)-0.1*zpix(ixc+2,iyc+1)-&
                                   0.1*zpix(ixc-2,iyc-1)-0.2*zpix(ixc-2,iyc)-0.1*zpix(ixc-2,iyc+1)-&
                                   0.1*zpix(ixc-1,iyc+2)-0.2*zpix(ixc,iyc+2)-0.1*zpix(ixc+1,iyc+2)-&
                                   0.1*zpix(ixc-1,iyc-2)-0.2*zpix(ixc,iyc-2)-0.1*zpix(ixc+1,iyc-2)&
                                    ))
                           endif
   
                           if (ikernel.eq.3) then
                              do i=-1,1
                                 do j=-1,1
                                    rd=abs(zpix(ix,iy)-zpix(ix+i,iy+j))
                                    if (rd.gt.l_diff_min) then
                                       rd=(rd-l_diff_min)/(l_diff_max-l_diff_min)
                                       l(ixl,iyl)=l(ixl,iyl)+min(rd,1.)
                                    endif
                                 enddo
                              enddo
                           endif
      
                           if (ikernel.eq.4) then
                              do i=-2,2
                                 do j=-2,2
                                    if (abs(i*j).ne.4) then
                                       rd=abs(zpix(ix,iy)-zpix(ix+i,iy+j))
                                       if (rd.gt.l_diff_min) then
                                          rd=abs(zpix(ix,iy)-zpix(ix+i,iy+j))
                                          rd=(rd-l_diff_min)/(l_diff_max-l_diff_min)
                                          l(ixl,iyl)=l(ixl,iyl)+min(rd,1.)
                                       endif
                                    endif
                                 enddo
                              enddo
                           endif
         
                           l(ixl,iyl)=min((l(ixl,iyl)-l_low)/(l_high-l_low),1.)
                           l(ixl,iyl)=max(l(ixl,iyl),0.)
                           if (l(ixl,iyl).gt.0.) rl=rl+1.
                           l_opacity_ave=l_opacity_ave+l(ixl,iyl)
                        enddo
                     enddo
                     if (rl.ge.6.) then
                        l_opacity=l_opacity_ave/6.
                     else
                        l_opacity=l(0,0)
                     endif
                     l_opacity=min(l_opacity,1.)
                     l_opacity=max(l_opacity,0.)
                  endif
                  ! ---- combine subunit outlines and derivative outlines ----
                  l_opacity=max(l_opacity,g_opacity)
               endif
   
               ! ***** CALCULATE THE TOTAL PIXEL INTENSITY *****
               ropacity=0.
               do icolor=1,3
                  rcolor= pfh*(pconetot*(colortype(type(atom(ix,iy)),icolor)))+(1.-pfh)*rfog(icolor) 
                  pix(ix,iy,icolor)=(1.-l_opacity)*rcolor
      
                  ! ----calculate pixel opacity
                  if (type(atom(ix,iy)).ne.0) ropacity=1.
                  pix(ix,iy,4)=max(ropacity,l_opacity)
      
               enddo
            enddo
   
            ! ***** output of a scan line *****
            ! ----- PPM format -----
            iscan=0
            do iout=1,iysize
               do ic=1,3
                  iscan=iscan+1
                  scanline(iscan)=int(pix(ix,iout,ic)*255.)
                  scanline(iscan)=min(scanline(iscan),255)
                  scanline(iscan)=max(scanline(iscan),0)
               enddo
            enddo
            write(8,'(20i4)') (scanline(if),if=1,iysize*3)
            ! -- write opacity
            !iscan=0
            !do iout=1,iysize
            !do ic=1,3
            !iscan=iscan+1
            !scanline(iscan)=int(pix(ix,iout,4)*255.)
            !scanline(iscan)=min(scanline(iscan),255)
            !scanline(iscan)=max(scanline(iscan),0)
            !enddo
            !enddo
            !write(9,'(20i4)') (scanline(if),if=1,iysize*3)
            ! ----- diagnostic ------
            if (int(ix/20)*20.eq.int((float(ix)/20.)*20.)) then
               write(6,'(65i1)') (int(pix(ix,iyo,1)*9.),iyo=1,iysize,20)
            endif
   
         enddo ! ix
   
      end select

   enddo

contains
!--------------------------------------------------------------------
   subroutine catenate(m1,m2)
      real*4 m1(4,4),m2(4,4),m(4,4)
      m(1,1)=m1(1,1)*m2(1,1)+m1(2,1)*m2(1,2)+m1(3,1)*m2(1,3)+m1(4,1)*m2(1,4)
      m(1,2)=m1(1,2)*m2(1,1)+m1(2,2)*m2(1,2)+m1(3,2)*m2(1,3)+m1(4,2)*m2(1,4)
      m(1,3)=m1(1,3)*m2(1,1)+m1(2,3)*m2(1,2)+m1(3,3)*m2(1,3)+m1(4,3)*m2(1,4)
      m(1,4)=m1(1,4)*m2(1,1)+m1(2,4)*m2(1,2)+m1(3,4)*m2(1,3)+m1(4,4)*m2(1,4)
      m(2,1)=m1(1,1)*m2(2,1)+m1(2,1)*m2(2,2)+m1(3,1)*m2(2,3)+m1(4,1)*m2(2,4)
      m(2,2)=m1(1,2)*m2(2,1)+m1(2,2)*m2(2,2)+m1(3,2)*m2(2,3)+m1(4,2)*m2(2,4)
      m(2,3)=m1(1,3)*m2(2,1)+m1(2,3)*m2(2,2)+m1(3,3)*m2(2,3)+m1(4,3)*m2(2,4)
      m(2,4)=m1(1,4)*m2(2,1)+m1(2,4)*m2(2,2)+m1(3,4)*m2(2,3)+m1(4,4)*m2(2,4)
      m(3,1)=m1(1,1)*m2(3,1)+m1(2,1)*m2(3,2)+m1(3,1)*m2(3,3)+m1(4,1)*m2(3,4)
      m(3,2)=m1(1,2)*m2(3,1)+m1(2,2)*m2(3,2)+m1(3,2)*m2(3,3)+m1(4,2)*m2(3,4)
      m(3,3)=m1(1,3)*m2(3,1)+m1(2,3)*m2(3,2)+m1(3,3)*m2(3,3)+m1(4,3)*m2(3,4)
      m(3,4)=m1(1,4)*m2(3,1)+m1(2,4)*m2(3,2)+m1(3,4)*m2(3,3)+m1(4,4)*m2(3,4)
      m(4,1)=m1(1,1)*m2(4,1)+m1(2,1)*m2(4,2)+m1(3,1)*m2(4,3)+m1(4,1)*m2(4,4)
      m(4,2)=m1(1,2)*m2(4,1)+m1(2,2)*m2(4,2)+m1(3,2)*m2(4,3)+m1(4,2)*m2(4,4)
      m(4,3)=m1(1,3)*m2(4,1)+m1(2,3)*m2(4,2)+m1(3,3)*m2(4,3)+m1(4,3)*m2(4,4)
      m(4,4)=m1(1,4)*m2(4,1)+m1(2,4)*m2(4,2)+m1(3,4)*m2(4,3)+m1(4,4)*m2(4,4)
      do j=1,4
         do i=1,4
            m1(i,j)=m(i,j)
         enddo
      enddo
      return
   end subroutine catenate
!--------------------------------------------------------------------
   subroutine clearmatrix(m)
      real*4 m(4,4)
      do i=1,4
         do j=1,4
            x=0.
            if (i.eq.j) x=1.
            m(i,j)=x
         enddo
      enddo
      return
   end subroutine clearmatrix
end
