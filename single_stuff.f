      subroutine task_init(myid_r,numtasks_r)
      implicit none
      integer myid_r,numtasks_r
      include 'msg_stuff.h'
      integer info, i
      integer npe_i, npe_j
      character*16 my_name
c
c single processor
      myid = 1
      mytid = 0
      tids(1) = mytid
      numtasks = 1
c ihost = hostnm_(my_name)
c print *,myid,mytid,' running on ',my_name
      myid_r = myid
      numtasks_r = numtasks
      return
      end
      subroutine task_end()
      implicit none
      include 'msg_stuff.h'
      integer info
      return
      end
      subroutine snd_real(dest, value, size, tag, info)
      implicit none
      integer dest, size, tag, info
      double precision value(size)
      include 'msg_stuff.h'
c print *,myid,'snd_real: ',dest,size,tag
      info = 0
      return
      end
      subroutine rcv_real(orig, value, size, tag, info)
      implicit none
      integer orig, size, tag, info
      double precision value(size)
      include 'msg_stuff.h'
      integer orig_tid
      info = 0
c print *,myid,'rcv_real: ',orig,size,tag,orig_tid
      return
      end
      subroutine snd_int(dest, value, size, tag, info)
      implicit none
      integer dest, size, tag, info
      integer value(size)
      include 'msg_stuff.h'
c print *,myid,'snd_int: ',dest,size,tag,value(1),tids(dest)
      info = 0
      return
      end
      subroutine rcv_int(orig, value, size, tag, info)
      implicit none
      integer orig, size, tag, info
      integer value(size)
      include 'msg_stuff.h'
      integer orig_tid
      info = 0
c print *,myid,'rcv_int: ',orig,size,tag,value(1),orig_tid
      return
      end
      subroutine bcast_real(value, size, tag, root, info)
      implicit none
      integer size, tag, root, info
      double precision value(size)
      include 'msg_stuff.h'
      if (numtasks .eq. 1) return
c print *,myid,'bcast_real: ',size,tag,root
      info = 0
      return
      end
      subroutine bcast_int(value, size, tag, root, info)
      implicit none
      integer size, tag, root, info
      integer value(size)
      include 'msg_stuff.h'
      if (numtasks .eq. 1) return
c print *,myid,'bcast_int: ',size,tag,root
      info = 0
      return
      end
      subroutine global_real_max(value)
      implicit none
      double precision value
      include 'msg_stuff.h'
      integer msgtag, i, info
      double precision y
      info = 0
      msgtag = 333
      if (myid .eq. 1) then
c receives values from nodes and find the max
         do i = 2, numtasks
            call rcv_real(-1, y, 1, msgtag, info)
            if (info.ne.0) print *,myid,'global_real_max: recv trouble'
            value = max(value,y)
         end do
      else
c tasks send their local max value to host:
         call snd_real(1, value, 1, msgtag, info)
         if (info.ne.0) print *,myid,'global_real_max: send trouble'
      endif
c master broadcasts global max
      info = 0
      msgtag = 334
      call bcast_real(value, 1, msgtag, 1, info)
      if (info.ne.0) print *,myid,'global_real_max: bcast trouble'
      return
      end
      subroutine global_real_sum(value)
      implicit none
      double precision value
      include 'msg_stuff.h'
      integer msgtag, i, info
      double precision y
      info = 0
      msgtag = 444
      if (myid .eq. 1) then
c receive values from nodes and sum them IN ORDER
         do i = 2, numtasks
            call rcv_real(i, y, 1, msgtag, info)
            if (info.ne.0) print *,myid,'global_real_sum: recv trouble'
            value = value + y
         end do
      else
c tasks send their local value to host:
         call snd_real(1, value, 1, msgtag, info)
         if (info.ne.0) print *,myid,'global_real_sum: send trouble'
      endif
c master broadcasts global sum
      info = 0
      msgtag = 445
      call bcast_real(value, 1, msgtag, 1, info)
      if (info.ne.0) print *,myid,'global_real_sum: bcast trouble'
      return
      end
      subroutine global_int_sum(value)
      implicit none
      integer value
      include 'msg_stuff.h'
      integer msgtag, i, info
      integer y
      info = 0
      msgtag = 444
      if (myid .eq. 1) then
c receive values from nodes and sum them IN ANY ORDER
         do i = 2, numtasks
            call rcv_real(-1, y, 1, msgtag, info)
            if (info.ne.0) print *,myid,'global_int_sum: recv trouble'
            value = value + y
         end do
      else
c tasks send their local value to host:
         call snd_real(1, value, 1, msgtag, info)
         if (info.ne.0) print *,myid,'global_int_sum: send trouble'
      endif
c master broadcasts global sum
      info = 0
      msgtag = 445
      call bcast_int(value, 1, msgtag, 1, info)
      if (info.ne.0) print *,myid,'global_int_sum: bcast trouble'
      return
      end
      subroutine barrier_sync()
      implicit none
      include 'msg_stuff.h'
      integer info
      if (numtasks .eq. 1) return
c print *,myid,'barrier_sync:'
      return
      end
