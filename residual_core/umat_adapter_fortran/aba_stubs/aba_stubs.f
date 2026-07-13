c=======================================================================
c  aba_stubs.f  --  standalone shims for the Abaqus solver utilities that
c  the Grilli UMAT references but that live inside Abaqus, not in the
c  (MIT) UMAT sources.
c
c  These are only ever CALLED at run time when the UMAT's optional
c  features are ON.  For the C3D8 crystal-plasticity replay we run with
c        twinon = 0   gndon = 0
c  so:
c    * MutexInit/Lock/Unlock  -> serial single-thread replay: no-ops.
c    * SMAIntArrayCreate/Access-> only used on the twin path: return 0
c      (never dereferenced because twinon=0).
c    * XIT                     -> Abaqus "abort job"; map to STOP.
c
c  They still have to be LINKED because the symbols are referenced in the
c  compiled object.  Signatures match SMAAspUserUtilities.hdr /
c  SMAAspUserArrays.hdr (integer(kind=4) ids, integer(kind=8) SMA
c  returns).  Under a real Abaqus/ifort build DO NOT link this file:
c  Abaqus provides these symbols itself.
c
c  Residual_Assembler, MIT.  Original work; no UMAT source copied.
c=======================================================================

c ---- Mutex (thread lock) no-ops --------------------------------------
      subroutine MutexInit(id)
      integer id
      return
      end

      subroutine MutexLock(id)
      integer id
      return
      end

      subroutine MutexUnlock(id)
      integer id
      return
      end

c ---- SMA global allocatable integer arrays (twin path only) ----------
c  Return a null handle.  With twinon=0 the returned pointer is never
c  associated with a Cray pointee that is read/written, so 0 is safe.
      integer*8 function SMAIntArrayCreate(id, isize, initval)
      integer id, isize, initval
      SMAIntArrayCreate = 0_8
      return
      end

      integer*8 function SMAIntArrayAccess(id)
      integer id
      SMAIntArrayAccess = 0_8
      return
      end

c ---- Abaqus "exit / abort analysis" ----------------------------------
      subroutine XIT
      write(*,'(A)') '*** UMAT called XIT (Abaqus abort). Stopping. ***'
      stop 1
      end
