      SUBROUTINE sweep (it,jt, kt, nm, isct, mm,mmo,mmi, mk, myid,
     1  hi, hj, hk, di, dj, dk, Phi, Phii,
     2  Src, Flux, Sigt,
     3  w,mu,eta,tsi, wmu,weta,wtsi, pn,
     4  north, south, east, west,
     5  Phiib,Phijb,Phikb,mdiag, Leakage,
     6  do_dsa,Face,it_dsa,jt_dsa,kt_dsa, do_fixup,nfixed,
     7  ibc,jbc,kbc,
     8  Phiibc,jt_ibc,kt_ibc,mm_ibc,
     9  Phijbc,it_jbc,kt_jbc,mm_jbc,
     A  Phikbc,it_kbc,jt_kbc,mm_kbc)


      IMPLICIT NONE
      INTEGER :: it, jt, kt, nm, isct, mm, mmo, mmi, mk, myid
      DOUBLE PRECISION :: hi(it), hj(jt), hk(kt)
      DOUBLE PRECISION :: di(it), dj(jt), dk(kt)
      DOUBLE PRECISION, DIMENSION(it) :: Phi
      DOUBLE PRECISION, DIMENSION(0:(it + 1)) :: Phii
      DOUBLE PRECISION, DIMENSION(it,jt,kt,nm) :: Src
      DOUBLE PRECISION, DIMENSION(it,jt,kt,nm) :: Flux
      DOUBLE PRECISION, DIMENSION(it,jt,kt) :: Sigt
      DOUBLE PRECISION :: w(mm), mu(mm), eta(mm), tsi(mm)
      DOUBLE PRECISION :: wmu(mm), weta(mm), wtsi(mm)
      DOUBLE PRECISION, DIMENSION(mm,nm,8) :: pn
      INTEGER :: north, south, east, west
      DOUBLE PRECISION, DIMENSION(jt,mk,mmi) :: Phiib
      DOUBLE PRECISION, DIMENSION(it,mk,mmi) :: Phijb
      DOUBLE PRECISION, DIMENSION(it,jt,mmi) :: Phikb
      INTEGER, DIMENSION(mmi) :: mdiag
      DOUBLE PRECISION, DIMENSION(6) :: Leakage
      LOGICAL :: do_dsa
      INTEGER :: it_dsa, jt_dsa, kt_dsa
      DOUBLE PRECISION, DIMENSION(it_dsa,jt_dsa,kt_dsa,3) :: Face
      LOGICAL :: do_fixup
      INTEGER :: nfixed
      INTEGER :: ibc, jbc, kbc
      INTEGER :: jt_ibc, kt_ibc, mm_ibc
      INTEGER :: it_jbc, kt_jbc, mm_jbc
      INTEGER :: it_kbc, jt_kbc, mm_kbc

      double precision phiibc(jt_ibc,kt_ibc,mm_ibc,0:1,0:1)
      double precision phijbc(it_jbc,kt_jbc,mm_jbc,0:1)
      double precision phikbc(it_kbc,jt_kbc,mm_kbc)
      INTEGER :: i, j, k, m, n, iq
      INTEGER :: idiag, ndiag, jkm, jk, mo, mi, mio
      INTEGER :: i0, i1, i2, i3, j0, j1, j2, j3, k0, k1, k2, k3
      INTEGER :: kk, kb, nk, lk
      INTEGER :: nib, njb
      INTEGER :: ns_snd, ns_rcv, ns_tag, ew_snd, ew_rcv, ew_tag, info
      DOUBLE PRECISION :: ci, cj, ck
      DOUBLE PRECISION :: phiir
      DOUBLE PRECISION :: leak
      DOUBLE PRECISION :: ql, dl, ti, tj, tk
      INTEGER :: ifixed, jfixed
      INTEGER :: nmess, mess
      COMMON / h / nmess,mess
C SWEEP FLOW
C
C DO iq=1,8                 ! octants
C
C  DO mo=1,mmo              ! angle pipelining loop
C
C   DO kk=1,kb              ! k-plane pipelining loop
C
C    RECV E/W               ! recv block I-inflows
C    RECV N/S               ! recv block J-inflows
C
C    DO idiag=1,jt+nk-1+mmi-1  ! JK-diagonals with MMI pipelining (block of work)
C     DO jkm=1,ndiag        ! I-lines on this diagonal
C
C      DO i=1,it            ! source (from Pn moments)
C      ENDDO
C
C      IF .NOT.do_fixups
C       DO i=i0,i1,i2       ! Sn eqn
C       ENDDO
C      ELSE
C       DO i=i0,i1,i2       ! Sn eqn w/ fixups
C       ENDDO
C      ENDIF
C
C      DO i=1,it            ! flux (Pn moments)
C      ENDDO
C
C      DO i=1,it            ! DSA face currents
C      ENDDO
C
C     ENDDO
C    ENDDO
C
C    SEND E/W               ! send block I-outflows
C    SEND N/S               ! send block J-outflows
C
C   ENDDO
C
C  ENDDO
C
C ENDDO
c     nmess = 0
c     mess = 0
      CALL indigo__record_f(2,108,Leakage(1),0,storage_size(Leakage))
      Leakage(1) = 0.0d+0
      CALL indigo__record_f(2,109,Leakage(2),0,storage_size(Leakage))
      Leakage(2) = 0.0d+0
      CALL indigo__record_f(2,110,Leakage(3),0,storage_size(Leakage))
      Leakage(3) = 0.0d+0
      CALL indigo__record_f(2,111,Leakage(4),0,storage_size(Leakage))
      Leakage(4) = 0.0d+0
      CALL indigo__record_f(2,112,Leakage(5),0,storage_size(Leakage))
      Leakage(5) = 0.0d+0
      CALL indigo__record_f(2,113,Leakage(6),0,storage_size(Leakage))
      Leakage(6) = 0.0d+0
      IF (do_dsa) THEN
      DO n = 1, 3
      DO k = 1, kt_dsa
      DO j = 1, jt_dsa
      DO i = 1, it_dsa
      CALL indigo__record_f(2,120,Face(i,j,k,n),1,storage_size(Face))
      Face(i,j,k,n) = 0.0d+0
      END DO
      END DO
      END DO
      END DO
      END IF
      nfixed = 0
c octant iq-loop
c
      DO iq = 1, 8
      CALL octant(i0,i1,i2,i3,j0,j1,j2,j3,k0,k1,k2,k3,iq,it,jt,kt)
      IF (i2 > 0) THEN
c  if i2 is > 0 (right sweep) recv phiib from west & send it east
      ew_rcv = west
      ew_snd = east
      ew_tag = 1000 * iq
      ELSE
c  if i2 is < 0 (left sweep)  recv phiib from east & send it west
      ew_rcv = east
      ew_snd = west
      ew_tag = 3000 * iq
      END IF
      IF (j2 > 0) THEN
c  if j2 is > 0 (up sweep)    recv phijb from south & send it north
      ns_rcv = south
      ns_snd = north
      ns_tag = 2000 * iq
      ELSE
c  if j2 is < 0  (down sweep) recv phijb from north & send it south
      ns_rcv = north
      ns_snd = south
      ns_tag = 4000 * iq
      END IF
      DO m = 1, mm
      CALL indigo__record_f(2,161,wmu(m),2,storage_size(wmu))
      CALL indigo__record_f(1,161,mu(m),3,storage_size(mu))
      CALL indigo__record_f(1,161,w(m),4,storage_size(w))
      wmu(m) = i2 * mu(m) * w(m)
      CALL indigo__record_f(2,162,weta(m),5,storage_size(weta))
      CALL indigo__record_f(1,162,eta(m),6,storage_size(eta))
      CALL indigo__record_f(1,162,w(m),4,storage_size(w))
      weta(m) = j2 * eta(m) * w(m)
      CALL indigo__record_f(2,163,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,163,tsi(m),8,storage_size(tsi))
      CALL indigo__record_f(1,163,w(m),4,storage_size(w))
      wtsi(m) = k2 * tsi(m) * w(m)
      END DO
c angle pipelining loop (batches of mmi angles)
c
      DO mo = 1, mmo
      mio = (mo - 1) * mmi
c K-inflows (k=k0 boundary)
c
      IF (k2 < 0 .OR. kbc .EQ. 0) THEN
      DO mi = 1, mmi
      DO j = 1, jt
      DO i = 1, it
      CALL indigo__record_f(2,177,Phikb(i,j,mi),9,storage_size(Phikb))
      Phikb(i,j,mi) = 0.0d+0
      END DO
      END DO
      END DO
      ELSE
      IF (do_dsa) THEN
      leak = 0.0
      k = k0 - k2
      DO mi = 1, mmi
      m = mi + mio
      DO j = 1, jt
      DO i = 1, it
      CALL indigo__record_f(2,189,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,189,Phikbc(i,j,m),10,storage_size(Phikbc))
      Phikb(i,j,mi) = Phikbc(i,j,m)
      CALL indigo__record_f(1,190,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,190,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,190,di(i),11,storage_size(di))
      CALL indigo__record_f(1,190,dj(j),12,storage_size(dj))
      leak = leak + wtsi(m) * Phikb(i,j,mi) * di(i) * dj(j)
      CALL indigo__record_f(2,192,Face(i,j,k + k3,3),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,192,Face(i,j,k + k3,3),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,192,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,192,Phikb(i,j,mi),9,storage_size(Phikb))
      Face(i,j,k + k3,3) = Face(i,j,k + k3,3) + wtsi(m) * Phikb(i,j,mi)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,197,Leakage(5),0,storage_size(Leakage))
      CALL indigo__record_f(1,197,Leakage(5),0,storage_size(Leakage))
      Leakage(5) = Leakage(5) + leak
      ELSE
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO j = 1, jt
      DO i = 1, it
      CALL indigo__record_f(2,204,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,204,Phikbc(i,j,m),10,storage_size(Phikbc))
      Phikb(i,j,mi) = Phikbc(i,j,m)
      CALL indigo__record_f(1,205,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,205,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,205,di(i),11,storage_size(di))
      CALL indigo__record_f(1,205,dj(j),12,storage_size(dj))
      leak = leak + wtsi(m) * Phikb(i,j,mi) * di(i) * dj(j)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,210,Leakage(5),0,storage_size(Leakage))
      CALL indigo__record_f(1,210,Leakage(5),0,storage_size(Leakage))
      Leakage(5) = Leakage(5) + leak
      END IF
      END IF
c k-plane pipelining loop (batches of mk-planes)
c
      kb = (kt + mk - 1) / mk
      DO kk = 1, kb
      IF (k2 > 0) THEN
      k0 = 1 + (kk - 1) * mk
      k1 = min(k0 + mk - 1,kt)
      nk = k1 - k0 + 1
      ELSE
      k0 = kt - (kk - 1) * mk
      k1 = max(k0 - mk + 1,1)
      nk = k0 - k1 + 1
      END IF
! this could be *nk* instead of *mk* if all phi{i,j}{b,bc}
! were dimensioned with mmi as the second dimension
      nib = jt * mk * mmi
      njb = it * mk * mmi
c I-inflows for block (i=i0 boundary)
c
      IF (ew_rcv .NE. 0) THEN
      CALL rcv_real(ew_rcv,Phiib,nib,ew_tag,info)
      ELSE
      IF (i2 < 0 .OR. ibc .EQ. 0) THEN
      DO mi = 1, mmi
      DO lk = 1, nk
      DO j = 1, jt
      CALL indigo__record_f(2,243,Phiib(j,lk,mi),13,storage_size(Phiib))
      Phiib(j,lk,mi) = 0.0d+0
      END DO
      END DO
      END DO
      ELSE
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO j = 1, jt
      CALL indigo__record_f(2,254,Phiib(j,lk,mi),13,storage_size(Phiib))
      CALL indigo__record_f(1,254,Phiibc(j,k,m,k3,j3),14,
     1 storage_size(Phiibc))
      Phiib(j,lk,mi) = Phiibc(j,k,m,k3,j3)
      CALL indigo__record_f(1,255,wmu(m),2,storage_size(wmu))
      CALL indigo__record_f(1,255,Phiib(j,lk,mi),13,storage_size(Phiib))
      CALL indigo__record_f(1,255,dj(j),12,storage_size(dj))
      CALL indigo__record_f(1,255,dk(k),15,storage_size(dk))
      leak = leak + wmu(m) * Phiib(j,lk,mi) * dj(j) * dk(k)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,260,Leakage(1),0,storage_size(Leakage))
      CALL indigo__record_f(1,260,Leakage(1),0,storage_size(Leakage))
      Leakage(1) = Leakage(1) + leak
      END IF
      END IF
      IF (do_dsa) THEN
      i = i0 - i2
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO j = 1, jt
      CALL indigo__record_f(2,270,Face(i + i3,j,k,1),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,270,Face(i + i3,j,k,1),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,270,wmu(m),2,storage_size(wmu))
      CALL indigo__record_f(1,270,Phiib(j,lk,mi),13,storage_size(Phiib))
      Face(i + i3,j,k,1) = Face(i + i3,j,k,1) + wmu(m) * Phiib(j,lk,mi)
      END DO
      END DO
      END DO
      END IF
c J-inflows for block (j=j0 boundary)
c
      IF (ns_rcv .NE. 0) THEN
      CALL rcv_real(ns_rcv,Phijb,njb,ns_tag,info)
      ELSE
      IF (j2 < 0 .OR. jbc .EQ. 0) THEN
      DO mi = 1, mmi
      DO lk = 1, nk
      DO i = 1, it
      CALL indigo__record_f(2,286,Phijb(i,lk,mi),16,storage_size(Phijb))
      Phijb(i,lk,mi) = 0.0d+0
      END DO
      END DO
      END DO
      ELSE
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO i = 1, it
      CALL indigo__record_f(2,297,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,297,Phijbc(i,k,m,k3),17,
     1 storage_size(Phijbc))
      Phijb(i,lk,mi) = Phijbc(i,k,m,k3)
      CALL indigo__record_f(1,298,weta(m),5,storage_size(weta))
      CALL indigo__record_f(1,298,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,298,di(i),11,storage_size(di))
      CALL indigo__record_f(1,298,dk(k),15,storage_size(dk))
      leak = leak + weta(m) * Phijb(i,lk,mi) * di(i) * dk(k)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,303,Leakage(3),0,storage_size(Leakage))
      CALL indigo__record_f(1,303,Leakage(3),0,storage_size(Leakage))
      Leakage(3) = Leakage(3) + leak
      END IF
      END IF
      IF (do_dsa) THEN
      j = j0 - j2
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO i = 1, it
      CALL indigo__record_f(2,313,Face(i,j + j3,k,2),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,313,Face(i,j + j3,k,2),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,313,weta(m),5,storage_size(weta))
      CALL indigo__record_f(1,313,Phijb(i,lk,mi),16,storage_size(Phijb))
      Face(i,j + j3,k,2) = Face(i,j + j3,k,2) + weta(m) * Phijb(i,lk,mi)
      END DO
      END DO
      END DO
      END IF
c JK-diagonals with MMI pipelined angles (block of work)
c
      DO mi = 1, mmi
      CALL indigo__record_f(2,323,mdiag(mi),18,storage_size(mdiag))
      mdiag(mi) = 0
      END DO
      DO idiag = 1, jt + nk - 1 + mmi - 1
      ndiag = 0
      DO mi = mmi, 2, -1
      CALL indigo__record_f(2,329,mdiag(mi),18,storage_size(mdiag))
      CALL indigo__record_f(1,329,mdiag(mi - 1),18,storage_size(mdiag))
      mdiag(mi) = mdiag(mi - 1)
      CALL indigo__record_f(1,330,mdiag(mi),18,storage_size(mdiag))
      ndiag = ndiag + mdiag(mi)
      END DO
      CALL indigo__record_f(2,332,mdiag(1),18,storage_size(mdiag))
      mdiag(1) = max(min(idiag,jt,nk,jt + nk - idiag),0)
      CALL indigo__record_f(1,333,mdiag(1),18,storage_size(mdiag))
      ndiag = ndiag + mdiag(1)
c I-lines on this diagonal (j/k/m triplets)
c
! DO PARALLEL
!   PRIVATE jkm,jk,mi,m,i,j,k,lk,n
!   PRIVATE ci,cj,ck
!   PRIVATE phi,phii,phiir
!   PRIVATE ql,dl,t,ti,tj,tk,ifixed,jfixed
!   SHARED mmi,mio,idiag,ndiag,mdiag
!   SHARED i0,j0,k0,i1,i2,j2,k2,it,jt,kt
!   SHARED iq,nm,pn,i3,j3,k3,mm,
!   SHARED w,mu,eta,tsi,wmu,weta,wtsi
!   SHARED src,sigt,flux,face
!   SHARED phiib,phijb,phikb
!   SHARED hi,hj,hk,di,dj,dk
!   SHARED do_fixup,do_dsa
!   SHARED nfixed
!   GUARD nfixed
      jfixed = 0
      DO jkm = 1, ndiag
      jk = jkm
      DO mi = 1, mmi - 1
      CALL indigo__record_f(1,356,mdiag(mi),18,storage_size(mdiag))
      IF (jk <= mdiag(mi)) GOTO 100
      CALL indigo__record_f(1,357,mdiag(mi),18,storage_size(mdiag))
      jk = jk - mdiag(mi)
      END DO
      mi = mmi
  100 CONTINUE
      m = mi + mio
      j = j0 + sign(min(idiag - mi + 1,jt) - 1 - jk + 1,j2)
      k = k0 + sign(max(idiag - mi + 1 - jt,0) + jk - 1,k2)
      lk = abs(k - k0) + 1
!
! Following 3 loops can also be used in place of above idiag/jkm loops,
! but then only the mmi loop can be task parallel.
!        DO mi = 1, mmi
!           m = mi + mio
!        DO lk = 1, nk
!           k = k0 + sign(lk-1,k2)
!        DO j = j0, j1, j2
      CALL indigo__record_f(1,374,tsi(m),8,storage_size(tsi))
      CALL indigo__record_f(1,374,hk(k),19,storage_size(hk))
      ck = tsi(m) * hk(k)
      CALL indigo__record_f(1,375,eta(m),6,storage_size(eta))
      CALL indigo__record_f(1,375,hj(j),20,storage_size(hj))
      cj = eta(m) * hj(j)
c I-inflow for this I-line
c
      CALL indigo__record_f(1,379,Phiib(j,lk,mi),13,storage_size(Phiib))
      phiir = Phiib(j,lk,mi)
      i = i0 - i2
      CALL indigo__record_f(2,381,Phii(i),21,storage_size(Phii))
      Phii(i) = phiir
c compute source from Pn moments (I-line)
      DO i = 1, it
      CALL indigo__record_f(2,385,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,385,Src(i,j,k,1),23,storage_size(Src))
      Phi(i) = Src(i,j,k,1)
      END DO
      DO n = 2, nm
      DO i = 1, it
      CALL indigo__record_f(2,389,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,389,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,389,pn(m,n,iq),24,storage_size(pn))
      CALL indigo__record_f(1,389,Src(i,j,k,n),23,storage_size(Src))
      Phi(i) = Phi(i) + pn(m,n,iq) * Src(i,j,k,n)
      END DO
      END DO
      IF (.NOT.do_fixup) THEN
c I-line recursion: without flux fixup
c
      DO i = i0, i1, i2
      CALL indigo__record_f(1,398,mu(m),3,storage_size(mu))
      CALL indigo__record_f(1,398,hi(i),25,storage_size(hi))
      ci = mu(m) * hi(i)
c balance equation - recursion on phiir
      CALL indigo__record_f(1,400,Sigt(i,j,k),26,storage_size(Sigt))
      dl = (Sigt(i,j,k) + ci + cj + ck)
      dl = 1.0 / dl
      CALL indigo__record_f(1,402,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,402,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,402,Phikb(i,j,mi),9,storage_size(Phikb))
      ql = (Phi(i) + ci * phiir + cj * Phijb(i,lk,mi) + ck * Phikb(i,j,
     1 mi))
      CALL indigo__record_f(2,404,Phi(i),22,storage_size(Phi))
      Phi(i) = ql * dl
c auxiliary equations (diamond)
      CALL indigo__record_f(1,406,Phi(i),22,storage_size(Phi))
      phiir = 2.0d+0 * Phi(i) - phiir
      CALL indigo__record_f(2,407,Phii(i),21,storage_size(Phii))
      Phii(i) = phiir
      CALL indigo__record_f(2,408,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,408,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,408,Phijb(i,lk,mi),16,storage_size(Phijb))
      Phijb(i,lk,mi) = 2.0d+0 * Phi(i) - Phijb(i,lk,mi)
      CALL indigo__record_f(2,409,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,409,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,409,Phikb(i,j,mi),9,storage_size(Phikb))
      Phikb(i,j,mi) = 2.0d+0 * Phi(i) - Phikb(i,j,mi)
      END DO
      ELSE
c I-line recursion: with flux fixup
c
      DO i = i0, i1, i2
      CALL indigo__record_f(1,417,mu(m),3,storage_size(mu))
      CALL indigo__record_f(1,417,hi(i),25,storage_size(hi))
      ci = mu(m) * hi(i)
c balance equation - recursion on phiir
      CALL indigo__record_f(1,419,Sigt(i,j,k),26,storage_size(Sigt))
      dl = (Sigt(i,j,k) + ci + cj + ck)
      ti = 1.0 / dl
      CALL indigo__record_f(1,421,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,421,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,421,Phikb(i,j,mi),9,storage_size(Phikb))
      ql = (Phi(i) + ci * phiir + cj * Phijb(i,lk,mi) + ck * Phikb(i,j,
     1 mi))
      CALL indigo__record_f(2,423,Phi(i),22,storage_size(Phi))
      Phi(i) = ql * ti
c auxiliary equations (diamond)
      CALL indigo__record_f(1,425,Phi(i),22,storage_size(Phi))
      ti = 2.0d+0 * Phi(i) - phiir
      CALL indigo__record_f(1,426,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,426,Phijb(i,lk,mi),16,storage_size(Phijb))
      tj = 2.0d+0 * Phi(i) - Phijb(i,lk,mi)
      CALL indigo__record_f(1,427,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,427,Phikb(i,j,mi),9,storage_size(Phikb))
      tk = 2.0d+0 * Phi(i) - Phikb(i,j,mi)
c fixup i,j, & k if negative
      ifixed = 0
  111 CONTINUE
      IF (ti < 0.0d+0) THEN
      dl = dl - ci
      ti = 1.0 / dl
      ql = ql - 0.5d+0 * ci * phiir
      CALL indigo__record_f(2,435,Phi(i),22,storage_size(Phi))
      Phi(i) = ql * ti
      ti = 0.0d+0
      IF (tj .NE. 0.0d+0) THEN
      CALL indigo__record_f(1,437,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,437,Phijb(i,lk,mi),16,storage_size(Phijb))
      tj = 2.0d+0 * Phi(i) - Phijb(i,lk,mi)
      END IF
      IF (tk .NE. 0.0d+0) THEN
      CALL indigo__record_f(1,438,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,438,Phikb(i,j,mi),9,storage_size(Phikb))
      tk = 2.0d+0 * Phi(i) - Phikb(i,j,mi)
      END IF
      ifixed = 1
      END IF
      IF (tj < 0.0d+0) THEN
      dl = dl - cj
      tj = 1.0 / dl
      CALL indigo__record_f(1,444,Phijb(i,lk,mi),16,storage_size(Phijb))
      ql = ql - 0.5d+0 * cj * Phijb(i,lk,mi)
      CALL indigo__record_f(2,445,Phi(i),22,storage_size(Phi))
      Phi(i) = ql * tj
      tj = 0.0d+0
      IF (tk .NE. 0.) THEN
      CALL indigo__record_f(1,447,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,447,Phikb(i,j,mi),9,storage_size(Phikb))
      tk = 2.0d+0 * Phi(i) - Phikb(i,j,mi)
      END IF
      IF (ti .NE. 0.) THEN
      CALL indigo__record_f(1,448,Phi(i),22,storage_size(Phi))
      ti = 2.0d+0 * Phi(i) - phiir
      END IF
      ifixed = 1
      GOTO 111
      END IF
      IF (tk < 0.0d+0) THEN
      dl = dl - ck
      tk = 1.0 / dl
      CALL indigo__record_f(1,455,Phikb(i,j,mi),9,storage_size(Phikb))
      ql = ql - 0.5d+0 * ck * Phikb(i,j,mi)
      CALL indigo__record_f(2,456,Phi(i),22,storage_size(Phi))
      Phi(i) = ql * tk
      tk = 0.0d+0
      IF (ti .NE. 0.0d+0) THEN
      CALL indigo__record_f(1,458,Phi(i),22,storage_size(Phi))
      ti = 2.0d+0 * Phi(i) - phiir
      END IF
      IF (tj .NE. 0.0d+0) THEN
      CALL indigo__record_f(1,459,Phi(i),22,storage_size(Phi))
      CALL indigo__record_f(1,459,Phijb(i,lk,mi),16,storage_size(Phijb))
      tj = 2.0d+0 * Phi(i) - Phijb(i,lk,mi)
      END IF
      ifixed = 1
      GOTO 111
      END IF
c
      phiir = ti
      CALL indigo__record_f(2,465,Phii(i),21,storage_size(Phii))
      Phii(i) = phiir
      CALL indigo__record_f(2,466,Phijb(i,lk,mi),16,storage_size(Phijb))
      Phijb(i,lk,mi) = tj
      CALL indigo__record_f(2,467,Phikb(i,j,mi),9,storage_size(Phikb))
      Phikb(i,j,mi) = tk
      jfixed = jfixed + ifixed
      END DO
      END IF
c compute flux Pn moments (I-line)
      DO i = 1, it
      CALL indigo__record_f(2,475,Flux(i,j,k,1),27,storage_size(Flux))
      CALL indigo__record_f(1,475,Flux(i,j,k,1),27,storage_size(Flux))
      CALL indigo__record_f(1,475,w(m),4,storage_size(w))
      CALL indigo__record_f(1,475,Phi(i),22,storage_size(Phi))
      Flux(i,j,k,1) = Flux(i,j,k,1) + w(m) * Phi(i)
      END DO
      DO n = 2, nm
      DO i = 1, it
      CALL indigo__record_f(2,479,Flux(i,j,k,n),27,storage_size(Flux))
      CALL indigo__record_f(1,479,Flux(i,j,k,n),27,storage_size(Flux))
      CALL indigo__record_f(1,479,pn(m,n,iq),24,storage_size(pn))
      CALL indigo__record_f(1,479,w(m),4,storage_size(w))
      CALL indigo__record_f(1,479,Phi(i),22,storage_size(Phi))
      Flux(i,j,k,n) = Flux(i,j,k,n) + pn(m,n,iq) * w(m) * Phi(i)
      END DO
      END DO
c compute DSA face currents (I-line)
      IF (do_dsa) THEN
      DO i = 1, it
      CALL indigo__record_f(2,487,Face(i + i3,j,k,1),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,487,Face(i + i3,j,k,1),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,487,wmu(m),2,storage_size(wmu))
      CALL indigo__record_f(1,487,Phii(i),21,storage_size(Phii))
      Face(i + i3,j,k,1) = Face(i + i3,j,k,1) + wmu(m) * Phii(i)
      CALL indigo__record_f(2,489,Face(i,j + j3,k,2),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,489,Face(i,j + j3,k,2),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,489,weta(m),5,storage_size(weta))
      CALL indigo__record_f(1,489,Phijb(i,lk,mi),16,storage_size(Phijb))
      Face(i,j + j3,k,2) = Face(i,j + j3,k,2) + weta(m) * Phijb(i,lk,mi)
      CALL indigo__record_f(2,491,Face(i,j,k + k3,3),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,491,Face(i,j,k + k3,3),1,
     1 storage_size(Face))
      CALL indigo__record_f(1,491,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,491,Phikb(i,j,mi),9,storage_size(Phikb))
      Face(i,j,k + k3,3) = Face(i,j,k + k3,3) + wtsi(m) * Phikb(i,j,mi)
      END DO
      END IF
c I-outflow for this I-line
c
      CALL indigo__record_f(2,498,Phiib(j,lk,mi),13,storage_size(Phiib))
      Phiib(j,lk,mi) = phiir
c diagonal I-line loops ( j/k/m triplets)
c
      END DO
      nfixed = nfixed + jfixed
      END DO
!
!        enddo
!        enddo
!        enddo
c block I-outflows (i=i1 boundary)
c
      IF (ew_snd .NE. 0) THEN
      CALL snd_real(ew_snd,Phiib,nib,ew_tag,info)
      ELSE
c              nmess = nmess + 1
c              mess = mess + nib
      IF (i2 < 0 .AND. ibc .NE. 0) THEN
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO j = 1, jt
      CALL indigo__record_f(2,524,Phiibc(j,k,m,k3,j3),
     1 14,storage_size(Phiibc))
      CALL indigo__record_f(1,524,Phiib(j,lk,mi),13,storage_size(Phiib))
      Phiibc(j,k,m,k3,j3) = Phiib(j,lk,mi)
      CALL indigo__record_f(1,525,wmu(m),2,storage_size(wmu))
      CALL indigo__record_f(1,525,Phiib(j,lk,mi),13,storage_size(Phiib))
      CALL indigo__record_f(1,525,dj(j),12,storage_size(dj))
      CALL indigo__record_f(1,525,dk(k),15,storage_size(dk))
      leak = leak + wmu(m) * Phiib(j,lk,mi) * dj(j) * dk(k)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,530,Leakage(1 + i3),0,
     1 storage_size(Leakage))
      CALL indigo__record_f(1,530,Leakage(1 + i3),0,
     1 storage_size(Leakage))
      Leakage(1 + i3) = Leakage(1 + i3) + leak
      ELSE
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO j = 1, jt
      CALL indigo__record_f(1,538,wmu(m),2,storage_size(wmu))
      CALL indigo__record_f(1,538,Phiib(j,lk,mi),13,storage_size(Phiib))
      CALL indigo__record_f(1,538,dj(j),12,storage_size(dj))
      CALL indigo__record_f(1,538,dk(k),15,storage_size(dk))
      leak = leak + wmu(m) * Phiib(j,lk,mi) * dj(j) * dk(k)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,543,Leakage(1 + i3),0,
     1 storage_size(Leakage))
      CALL indigo__record_f(1,543,Leakage(1 + i3),0,
     1 storage_size(Leakage))
      Leakage(1 + i3) = Leakage(1 + i3) + leak
      END IF
      END IF
c block J-outflows (j=j1 boundary)
c
      IF (ns_snd .NE. 0) THEN
      CALL snd_real(ns_snd,Phijb,njb,ns_tag,info)
      ELSE
c              nmess = nmess + 1
c              mess = mess + njb
      IF (j2 < 0 .AND. jbc .NE. 0) THEN
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO i = 1, it
      CALL indigo__record_f(2,561,Phijbc(i,k,
     1 m,k3),17,storage_size(Phijbc))
      CALL indigo__record_f(1,561,Phijb(i,lk,mi),16,storage_size(Phijb))
      Phijbc(i,k,m,k3) = Phijb(i,lk,mi)
      CALL indigo__record_f(1,562,weta(m),5,storage_size(weta))
      CALL indigo__record_f(1,562,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,562,di(i),11,storage_size(di))
      CALL indigo__record_f(1,562,dk(k),15,storage_size(dk))
      leak = leak + weta(m) * Phijb(i,lk,mi) * di(i) * dk(k)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,567,Leakage(3 + j3),0,
     1 storage_size(Leakage))
      CALL indigo__record_f(1,567,Leakage(3 + j3),0,
     1 storage_size(Leakage))
      Leakage(3 + j3) = Leakage(3 + j3) + leak
      ELSE
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO lk = 1, nk
      k = k0 + sign(lk - 1,k2)
      DO i = 1, it
      CALL indigo__record_f(1,575,weta(m),5,storage_size(weta))
      CALL indigo__record_f(1,575,Phijb(i,lk,mi),16,storage_size(Phijb))
      CALL indigo__record_f(1,575,di(i),11,storage_size(di))
      CALL indigo__record_f(1,575,dk(k),15,storage_size(dk))
      leak = leak + weta(m) * Phijb(i,lk,mi) * di(i) * dk(k)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,580,Leakage(3 + j3),0,
     1 storage_size(Leakage))
      CALL indigo__record_f(1,580,Leakage(3 + j3),0,
     1 storage_size(Leakage))
      Leakage(3 + j3) = Leakage(3 + j3) + leak
      END IF
      END IF
c k-plane pipelining loop (batches of mk-planes)
c
      END DO
c K-outflows (k=k1 boundary)
c
      IF (k2 < 0 .AND. kbc .NE. 0) THEN
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO j = 1, jt
      DO i = 1, it
      CALL indigo__record_f(2,596,Phikbc(i,j,m),10,storage_size(Phikbc))
      CALL indigo__record_f(1,596,Phikb(i,j,mi),9,storage_size(Phikb))
      Phikbc(i,j,m) = Phikb(i,j,mi)
      CALL indigo__record_f(1,597,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,597,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,597,di(i),11,storage_size(di))
      CALL indigo__record_f(1,597,dj(j),12,storage_size(dj))
      leak = leak + wtsi(m) * Phikb(i,j,mi) * di(i) * dj(j)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,602,Leakage(5 + k3),0,
     1 storage_size(Leakage))
      CALL indigo__record_f(1,602,Leakage(5 + k3),0,
     1 storage_size(Leakage))
      Leakage(5 + k3) = Leakage(5 + k3) + leak
      ELSE
      leak = 0.0
      DO mi = 1, mmi
      m = mi + mio
      DO j = 1, jt
      DO i = 1, it
      CALL indigo__record_f(1,609,wtsi(m),7,storage_size(wtsi))
      CALL indigo__record_f(1,609,Phikb(i,j,mi),9,storage_size(Phikb))
      CALL indigo__record_f(1,609,di(i),11,storage_size(di))
      CALL indigo__record_f(1,609,dj(j),12,storage_size(dj))
      leak = leak + wtsi(m) * Phikb(i,j,mi) * di(i) * dj(j)
      END DO
      END DO
      END DO
      CALL indigo__record_f(2,614,Leakage(5 + k3),0,
     1 storage_size(Leakage))
      CALL indigo__record_f(1,614,Leakage(5 + k3),0,
     1 storage_size(Leakage))
      Leakage(5 + k3) = Leakage(5 + k3) + leak
      END IF
c angle pipelining loop (batches of mmi angles)
c
      END DO
c octant loop
c
      END DO
      RETURN
      END SUBROUTINE 

      SUBROUTINE indigo__create_map()
      CALL indigo__write_idx_f("Leakage",7)
      CALL indigo__write_idx_f("Face",4)
      CALL indigo__write_idx_f("wmu",3)
      CALL indigo__write_idx_f("mu",2)
      CALL indigo__write_idx_f("w",1)
      CALL indigo__write_idx_f("weta",4)
      CALL indigo__write_idx_f("eta",3)
      CALL indigo__write_idx_f("wtsi",4)
      CALL indigo__write_idx_f("tsi",3)
      CALL indigo__write_idx_f("Phikb",5)
      CALL indigo__write_idx_f("Phikbc",6)
      CALL indigo__write_idx_f("di",2)
      CALL indigo__write_idx_f("dj",2)
      CALL indigo__write_idx_f("Phiib",5)
      CALL indigo__write_idx_f("Phiibc",6)
      CALL indigo__write_idx_f("dk",2)
      CALL indigo__write_idx_f("Phijb",5)
      CALL indigo__write_idx_f("Phijbc",6)
      CALL indigo__write_idx_f("mdiag",5)
      CALL indigo__write_idx_f("hk",2)
      CALL indigo__write_idx_f("hj",2)
      CALL indigo__write_idx_f("Phii",4)
      CALL indigo__write_idx_f("Phi",3)
      CALL indigo__write_idx_f("Src",3)
      CALL indigo__write_idx_f("pn",2)
      CALL indigo__write_idx_f("hi",2)
      CALL indigo__write_idx_f("Sigt",4)
      CALL indigo__write_idx_f("Flux",4)
      END SUBROUTINE 

