      PROGRAM driver
c **********************************************************************
c *                                                                    *
c *  copyright, 1993, the regents of the university of california.     *
c *  this software was produced under a u. s. government contract      *
c *  (w-7405-eng-36) by the los alamos national laboratory, which is   *
c *  operated by the university of california for the u. s. department *
c *  of energy. the u. s. government is licensed to use, reproduce,    *
c *  and distribute this software. permission is granted to the public *
c *  to copy and use this software without charge, provided that this  *
c *  notice and any statement of authorship are reproduced on all      *
c *  copies. neither the government nor the university makes any       *
c *  warranty, express or implied, or assumes any liability            *
c *  responsibility for the use of this software.                      *
c *                                                                    *
c **********************************************************************
c
c  this program is an explict decomposition version of the 3d sweeper
c
      IMPLICIT NONE
      INTEGER :: it_g, jt_g, it, jt, kt, mm, mmo, mmi, isct
      DOUBLE PRECISION :: dx, dy, dz, epsi, err
      INTEGER :: mk, isn, nm
      INTEGER :: iprint
      LOGICAL :: do_dsa
      INTEGER :: it_dsa, jt_dsa, kt_dsa
      INTEGER :: ibc, jbc, kbc
      INTEGER :: jt_ibc, kt_ibc, mm_ibc
      INTEGER :: it_jbc, kt_jbc, mm_jbc
      INTEGER :: it_kbc, jt_kbc, mm_kbc
      INTEGER :: myid, numtasks, npe_i, npe_j, mype_i, mype_j
      INTEGER :: istart, jstart
      INTEGER :: north, south, east, west
      INTEGER :: its
      DOUBLE PRECISION :: t0, t1, et0, et1, grind
      INTEGER :: kb, nmsgs
      DOUBLE PRECISION :: mem, effic, effic1, effic2
      INTEGER :: i, mi, ndiag, icpu, ncpu
      INTEGER :: ifixups, nfixups
c
c  initialize the tasking system and spawn the daughter tasks:
c    the processor grid is npe_i X npe_j
c    send pvm initialization data to other tasks:
c    set position of each task in the processor grid (mype_i, mype_j):
c    determine the subgrid (it by jt):
c    determine the neighbors for each processor task:
c
      CALL indigo__init(1,0)
      CALL indigo__create_map()
      CALL task_init(myid,numtasks)
c
c read input on master node and broadcast to slaves
c
      CALL read_input(it_g,jt_g,kt,mm,isct,dx,dy,dz,epsi,iprint,do_dsa,ifixups,ibc,jbc,kbc,mk,isn,nm,myid,numtasks,npe_i,npe_j,mmo,mmi,ncpu)
c
c determine the decomposition
c
      CALL decomp(it_g,jt_g,it,jt,istart,jstart,kt,myid,numtasks,npe_i,npe_j,mype_i,mype_j,north,south,east,west)
c
c conditional array sizes
c
      IF (do_dsa) THEN
      it_dsa = it + 1
      jt_dsa = jt + 1
      kt_dsa = kt + 1
      ELSE
      it_dsa = 1
      jt_dsa = 1
      kt_dsa = 1
      END IF
      IF (ibc .NE. 0) THEN
      jt_ibc = jt
      kt_ibc = kt
      mm_ibc = mm
      ELSE
      jt_ibc = 1
      kt_ibc = 1
      mm_ibc = 1
      END IF
      IF (jbc .NE. 0) THEN
      it_jbc = it
      kt_jbc = kt
      mm_jbc = mm
      ELSE
      it_jbc = 1
      kt_jbc = 1
      mm_jbc = 1
      END IF
      IF (kbc .NE. 0) THEN
      it_kbc = it
      jt_kbc = jt
      mm_kbc = mm
      ELSE
      it_kbc = 1
      jt_kbc = 1
      mm_kbc = 1
      END IF
      IF (myid .EQ. 1) THEN
      PRINT *, 'SWEEP3D - Method 5 -',' Pipelined Wavefront with Line-Recursion'
      PRINT *, 'Version 2.2b'
      PRINT 100, isn,isct,mm,nm,it_g,jt_g,kt
  100 format ( ' S',i1,'P',i1,3x,'-',1x,i2,' angles/octant,',i3,' moments',/,' global grid: ',i3,'x',i3,'x',i3 )
      PRINT *, numtasks,'domains   - ',npe_i,'x',npe_j,'decomposition'
      mem = float(it) * float(jt) * float(kt) * (3 + (isct + 1) + 2 * nm)
      IF (do_dsa) THEN
      mem = mem + float(it_dsa) * float(jt_dsa) * float(kt_dsa) * 3
      END IF
      kb = (kt + mk - 1) / mk
      PRINT *, kb * mmo,'domain pipelined blocks -',mk,'k-planes by',mmi,'angles each'
      mem = mem * 8 / 1.0d+6
      PRINT 101, mem
  101 format ( ' estimated memory usage per domain:',f7.1,' MB' )
      IF (npe_i .EQ. 1 .AND. npe_j .EQ. 1) THEN
      nmsgs = 0
      ELSE IF (npe_i .EQ. 1) THEN
      nmsgs = 2 * (2 * mmo) * (1 + 1 + 0 + 0) + (npe_j - 2) * (2 * mmo) * (1 + 1 + 1 + 1)
      ELSE IF (npe_j .EQ. 1) THEN
      nmsgs = 2 * (2 * mmo) * (1 + 1 + 0 + 0) + (npe_i - 2) * (2 * mmo) * (1 + 1 + 1 + 1)
      ELSE 
      nmsgs = 4 * (2 * mmo) * (2 + 1 + 0 + 1) + 2 * (npe_i - 2) * (2 * mmo) * (2 + 1 + 1 + 2) + 2 * (npe_j - 2) * (2 * mmo) * (2 + 1 + 1 + 2) + (npe_i - 2) * (npe_j - 2) * (2 * mmo) * (2 + 2 + 2 + 2)
      END IF
      nmsgs = nmsgs * kb
      PRINT *, nmsgs,'global messages per iteration'
      effic1 = float(8 * mmo * kb) / float(8 * mmo * kb + 2 * (npe_i - 1) + 4 * (npe_j - 1))
      PRINT 102, 100.0 * effic1
  102 format ( f6.2,'% domain parallel efficiency -',' due to decomposition & blocking' )
      icpu = 0
      DO i = 1, jt + mk - 1 + mmi - 1
      ndiag = 0
      DO mi = 1, mmi
      ndiag = ndiag + max(min(i - mi + 1,jt,mk,jt + mk - i + mi - 1),0)
      END DO
      icpu = icpu + (ndiag + ncpu - 1) / ncpu
      END DO
      effic2 = float(mmi * jt * mk) / float(icpu * ncpu)
      PRINT 103, 100.0 * effic2,ncpu
  103 format ( f6.2,'% multitasking efficiency on',i4,' processors' )
      PRINT 104, 100.0 * effic1 * effic2,npe_i * npe_j,ncpu
  104 format ( f6.2,'% combined efficiency on a cluster of ',i4,i4,'-way SMPs' )
      IF (do_dsa) THEN
      PRINT *, 'DSA leakage calculation: ON'
      ELSE
      PRINT *, 'DSA leakage calculation: OFF'
      END IF
      IF (ifixups > 0) THEN
      PRINT *, 'Flux fixups: ON (always)'
      ELSE IF (ifixups .EQ. 0) THEN
      PRINT *, 'Flux fixups: OFF'
      ELSE 
      PRINT *, 'Flux fixups: ON (after',-ifixups,'iterations)'
      END IF
      END IF
      CALL inner_auto(it,jt,kt,nm,isct,mm,mmo,mmi,mk,dx,dy,dz,epsi,err,iprint,do_dsa,it_g,jt_g,istart,jstart,myid,numtasks,north,south,east,west,its,t0,t1,et0,et1,it_dsa,jt_dsa,kt_dsa,ifixups,nfixups,ibc,jbc,kbc,jt_ibc,kt_ibc,mm_ibc,it_jbc,kt_jbc,mm_jbc,it_kbc,jt_kbc,mm_kbc)
      IF (myid .EQ. 1) THEN
      PRINT *, 'CPU     time was: ',t1 - t0
      PRINT *, 'Elapsed time was: ',et1 - et0
      grind = 1.0e6 * (t1 - t0) / (real(its) * real(it_g) * real(jt_g) * real(kt) * real(mm) * real(8))
      PRINT '('' CPU grind time:  '',1pg10.3)', grind
      grind = 1.0e6 * (et1 - et0) / (real(its) * real(it_g) * real(jt_g) * real(kt) * real(mm) * real(8))
      PRINT '('' Wall grind time: '',1pg10.3)', grind
      END IF
      CALL task_end()
      STOP 
      END PROGRAM 

