/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
#ifndef __uSWSR_PTR_BUFFER_HPP_
#define __uSWSR_PTR_BUFFER_HPP_
/* ***************************************************************************
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License version 3 as 
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *
 ****************************************************************************
 */

/* 
 * Single-Writer/Single-Reader (SWSR) lock-free (wait-free) unbounded 
 * FIFO queue.  No lock is needed around pop and push methods!!
 * 
 * The key idea underneath the implementation is quite simple: 
 * the unbounded queue is based on a pool of wait-free 
 * SWSR circular buffers (see buffer.hpp).  The pool of buffers 
 * automatically grows and shrinks on demand. The implementation 
 * of the pool of buffers carefully try to minimize the impact of 
 * dynamic memory allocation/deallocation by using caching strategies. 
 *
 * More information about the uSWSR_Ptr_Buffer implementation and 
 * correctness proof can be found in the following report:
 *
 * Massimo Torquati, "Single-Producer/Single-Consumer Queue on Shared Cache 
 * Multi-Core Systems", TR-10-20, Computer Science Department, University
 * of Pisa Italy,2010
 * ( http://compass2.di.unipi.it/TR/Files/TR-10-20.pdf.gz )
 *
 *
 *
 * IMPORTANT:
 *
 * This implementation has been optimized for 1 producer and 1 consumer. 
 * If you need to use more producers and/or more consumers you have 
 * several possibilities (top-down order):
 *  1. to use an high level construct like the farm skeleton and compositions 
 *     of multiple farms (in other words, use FastFlow ;-) ).
 *  2. to use one of the implementations in the MPMCqueues.hpp file
 *  3. to use the SWSR_Ptr_Buffer but with the mp_push and the mc_pop methods 
 *     both protected by (spin-)locks in order to protect the internal data 
 *     structures.
 *
 *  
 */
#include <assert.h>
#include <ff/dynqueue.hpp>
#include <ff/buffer.hpp>
#include <ff/spin-lock.hpp>
#if defined(HAVE_ATOMIC_H)
#include <asm/atomic.h>
#else
#include <ff/atomic/atomic.h>
#endif
namespace ff {

/* Do not change the following define unless you know what you're doing */
#define INTERNAL_BUFFER_T SWSR_Ptr_Buffer  /* Lamport_Buffer */

class BufferPool {
public:
    BufferPool(int cachesize)
        :inuse(cachesize),bufcache(cachesize) {
        bufcache.init();
#if defined(UBUFFER_STATS)
        miss=0;hit=0;
#endif
    }
    
    ~BufferPool() {
        union { INTERNAL_BUFFER_T * b1; void * b2;} p;

        while(inuse.pop(&p.b2)) {
            p.b1->~INTERNAL_BUFFER_T();
            free(p.b2);
        }
        while(bufcache.pop(&p.b2)) {
            p.b1->~INTERNAL_BUFFER_T();	    
            free(p.b2);
        }
    }
    
    inline INTERNAL_BUFFER_T * const next_w(unsigned long size)  { 
        union { INTERNAL_BUFFER_T * buf; void * buf2;} p;
        if (!bufcache.pop(&p.buf2)) {
#if defined(UBUFFER_STATS)
            ++miss;
#endif
            p.buf = (INTERNAL_BUFFER_T*)malloc(sizeof(INTERNAL_BUFFER_T));
            new (p.buf) INTERNAL_BUFFER_T(size);
#if defined(uSWSR_MULTIPUSH)        
            if (!p.buf->init(true)) return NULL;
#else
            if (!p.buf->init()) return NULL;
#endif
        }
#if defined(UBUFFER_STATS)
        else  ++hit;
#endif  
        inuse.push(p.buf);
        return p.buf;
    }
    
    inline INTERNAL_BUFFER_T * const next_r()  { 
        union { INTERNAL_BUFFER_T * buf; void * buf2;} p;
        return (inuse.pop(&p.buf2)? p.buf : NULL);
    }
    
    inline void release(INTERNAL_BUFFER_T * const buf) {
        buf->reset();
        if (!bufcache.push(buf)) {
            buf->~INTERNAL_BUFFER_T();	    
            free(buf);
        }
    }

#if defined(UBUFFER_STATS)
    inline unsigned long readPoolMiss() {
        unsigned long m = miss;
        miss = 0;
        return m;
    }
    inline unsigned long readPoolHit() {
        unsigned long h = hit;
        hit = 0;
        return h;
    }
#endif

private:
#if defined(UBUFFER_STATS)
    unsigned long      miss,hit;
    long padding1[longxCacheLine-2];    
#endif

    dynqueue           inuse;
    INTERNAL_BUFFER_T  bufcache;
};
    

// unbounded buffer based on the INTERNAL_BUFFER_T 
class uSWSR_Ptr_Buffer {
private:
    enum {CACHE_SIZE=32};

#if defined(uSWSR_MULTIPUSH)
    enum { MULTIPUSH_BUFFER_SIZE=16};

    inline bool multipush() {
        if (buf_w->multipush(multipush_buf,MULTIPUSH_BUFFER_SIZE)) {
            mcnt=0; 
            return true;
        }

        if (fixedsize) return false;
        // try to get a new buffer             
        INTERNAL_BUFFER_T * t = pool.next_w(size);
        assert(t); // if (!t) return false; // EWOULDBLOCK
        buf_w = t;
        buf_w->multipush(multipush_buf,MULTIPUSH_BUFFER_SIZE);
        mcnt=0;
#if defined(UBUFFER_STATS)
        atomic_long_inc(&numBuffers);
#endif
        return true;
    }
#endif

public:
    uSWSR_Ptr_Buffer(unsigned long n, const bool fixedsize=false):
        buf_r(0),buf_w(0),size(n),fixedsize(fixedsize),
        pool(CACHE_SIZE) {
        init_unlocked(P_lock); init_unlocked(C_lock);
#if defined(UBUFFER_STATS)
        atomic_long_set(&numBuffers,0);
#endif
    }
    
    ~uSWSR_Ptr_Buffer() {
        if (buf_r) {
            buf_r->~INTERNAL_BUFFER_T();
            free(buf_r);
        }
        // buf_w either is equal to buf_w or is freed by BufferPool destructor
    }
    
    bool init() {
        if (buf_w || buf_r) return false;
#if defined(uSWSR_MULTIPUSH)
        if (size<=MULTIPUSH_BUFFER_SIZE) return false;
#endif
        buf_r = (INTERNAL_BUFFER_T*)::malloc(sizeof(INTERNAL_BUFFER_T));
        if (!buf_r) return false;
        new ((void *)buf_r) INTERNAL_BUFFER_T(size);
#if defined(uSWSR_MULTIPUSH)        
        if (!buf_r->init(true)) return false;
#else
        if (!buf_r->init()) return false;
#endif
        buf_w = buf_r;
        return true;
    }
    
    /* return true if the buffer is empty */
    inline bool empty() {
        return buf_r->empty();
    }
    
    /* return true if there is at least one room in the buffer */
    inline bool available()   { 
        return buf_w->available();
    }
 
    /* push the input value into the queue. 
     * If fixedsize has been set, this method may
     * returns false, this means EWOULDBLOCK 
     * and the call should be retried.  
     */
    inline bool push(void * const data) {
        /* NULL values cannot be pushed in the queue */
        if (!data || !buf_w) return false;
        
        if (!available()) {

            if (fixedsize) return false;

            // try to get a new buffer             
            INTERNAL_BUFFER_T * t = pool.next_w(size);
            assert(t); //if (!t) return false; // EWOULDBLOCK
            buf_w = t;
#if defined(UBUFFER_STATS)
            atomic_long_inc(&numBuffers);
#endif
        }
        //DBG(assert(buf_w->push(data)); return true;);
        buf_w->push(data);
        return true;
    }

    /* multi-producers method. lock protected.
     *
     */
    inline bool mp_push(void *const data) {
        spin_lock(P_lock);
        bool r=push(data);
        spin_unlock(P_lock);
        return r;
    }

#if defined(uSWSR_MULTIPUSH)
    /* massimot: experimental code
     *
     * This method provides the same interface of the push one but 
     * uses the multipush method to provide a batch of items to
     * the consumer thus ensuring better cache locality and 
     * lowering the cache trashing.
     */
    inline bool mpush(void * const data) {
        if (!data) return false;
        
        if (mcnt==MULTIPUSH_BUFFER_SIZE) 
            return multipush();        

        multipush_buf[mcnt++]=data;

        if (mcnt==MULTIPUSH_BUFFER_SIZE) return multipush();

        return true;
    }

    inline bool flush() {
        return (mcnt ? multipush() : true);
    }
#endif /* uSWSR_MULTIPUSH */
    
    inline bool  pop(void ** data) {
        if (!data || !buf_r) return false;
        if (buf_r->empty()) { // current buffer is empty
            if (buf_r == buf_w) return false;
            if (buf_r->empty()) { // we have to check again
                INTERNAL_BUFFER_T * tmp = pool.next_r();
                if (tmp) {
                    // there is another buffer, release the current one 
                    pool.release(buf_r); 
                    buf_r = tmp;                    

#if defined(UBUFFER_STATS)
                    atomic_long_dec(&numBuffers);
#endif
                }
            }
        }
        //DBG(assert(buf_r->pop(data)); return true;);
        return buf_r->pop(data);
    }    

#if defined(UBUFFER_STATS)
    inline unsigned long queue_status() {
        return (unsigned long)atomic_long_read(&numBuffers);
    }
    inline unsigned long readMiss() {
        return pool.readPoolMiss();
    }
    inline unsigned long readHit() {
        return pool.readPoolHit();
    }
#endif

  /* multi-consumers method. lock protected.
     *
     */
    inline bool mc_pop(void ** data) {
        spin_lock(C_lock);
        bool r=pop(data);
        spin_unlock(C_lock);
        return r;
    }
    
    /* NOTE: this is not the real queue length
     * but just a rough estimation.
     *
     */
    inline unsigned long length() const {
        unsigned long len = buf_r->length();
        if (buf_r == buf_w) return len;
        return len+buf_w->length();
    }


private:
    // Padding is required to avoid false-sharing between 
    // core's private cache
    INTERNAL_BUFFER_T * buf_r;
    long padding1[longxCacheLine-1];
    INTERNAL_BUFFER_T * buf_w;
    long padding2[longxCacheLine-1];

    /* ----- two-lock used only in the mp_push and mp_pop methods ------- */
    lock_t P_lock;
    long padding3[longxCacheLine-sizeof(lock_t)];
    lock_t C_lock;
    long padding4[longxCacheLine-sizeof(lock_t)];
    /* -------------------------------------------------------------- */
#if defined(UBUFFER_STATS)
    atomic_long_t numBuffers;
#endif

#if defined(uSWSR_MULTIPUSH)
    /* massimot: experimental code (see multipush)
     *
     */
    // local multipush buffer used by the mpush method
    void  * multipush_buf[MULTIPUSH_BUFFER_SIZE];
    int     mcnt;
#endif
    const unsigned long	size;
    const bool			fixedsize;
    BufferPool			pool;
};
    
} // namespace


#endif /* __uSWSR_PTR_BUFFER_HPP_ */
