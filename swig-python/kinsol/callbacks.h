/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This header file defines 'director' classes which allow python
 * functions to be provided as the residual and Jacboian.
 * -----------------------------------------------------------------*/

#include <iostream>
#include <functional>

class KINSysPyFn {
  public:
    virtual ~KINSysPyFn() {}
    virtual int actual_sysfun(N_Vector y, N_Vector g, void *udata) { return -1; }
};

class KINSysFnCaller {
  private:
    KINSysPyFn *sysfn_;
    N_Vector y_;
    N_Vector g_;
    void *udata_;

  public:
    KINSysFnCaller(): sysfn_(nullptr), y_(nullptr), g_(nullptr), udata_(nullptr) {}
    ~KINSysFnCaller() { cleanup(); }
    
    void cleanup()
    { 
      delete sysfn_;
      sysfn_ = nullptr;
    }
    
    void setFn(KINSysPyFn *cb)
    {
      cleanup();
      sysfn_ = cb;
    }
    
    void setArgs(N_Vector y, N_Vector g, void *udata = nullptr)
    {
      y_ = y;
      g_ = g;
      udata_ = udata;
    }
    
    int call() 
    {
      if(sysfn_) return sysfn_->actual_sysfun(y_, g_, udata_);
      else       return -1;
    }
};

int KINPyInterfaceSysFn(N_Vector y, N_Vector g, void *user_data)
{
  KINSysFnCaller *caller = static_cast<KINSysFnCaller*>(user_data);
  caller->setArgs(y, g);
  return caller->call();
}

int KINInitPy(void *kmem, KINSysFnCaller *caller, N_Vector y)
{
  int flag;
  flag = KINInit(kmem, KINPyInterfaceSysFn, y);
  if (flag < 0) return flag;

  /* hijack user data so we can provide the caller */
  flag = KINSetUserData(kmem, static_cast<void*>(caller));
  return flag;
}
