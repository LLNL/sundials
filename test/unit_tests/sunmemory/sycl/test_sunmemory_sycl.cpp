/*------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------*/

#include <CL/sycl.hpp>
#include <iostream>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_types.h>
#include <sunmemory/sunmemory_sycl.h>

int test_instance(SUNMemoryHelper helper, SUNMemoryType mem_type,
                  bool print_test_status)
{
  // Create an in-order GPU queue
#if SYCL_LANGUAGE_VERSION >= 2020
  sycl::queue myQueue(sycl::gpu_selector_v,
                      sycl::property_list{sycl::property::queue::in_order{}});
#else
  sycl::gpu_selector selector;
  sycl::queue myQueue(selector,
                      sycl::property_list{sycl::property::queue::in_order{}});
#endif

  // Try and allocate some memory
  const int N                 = 8;
  const size_t bytes_to_alloc = N * sizeof(sunrealtype);
  SUNMemory some_memory       = nullptr;

  if (print_test_status) std::cout << "  SUNMemoryHelper_Alloc... \n";
  int retval = SUNMemoryHelper_Alloc(helper, &some_memory, bytes_to_alloc,
                                     mem_type, static_cast<void*>(&myQueue));
  if (retval)
  {
    if (print_test_status) std::cout << "  SUNMemoryHelper_Alloc... FAILED\n";
    return -1;
  }

  // Write to the memory
  sunrealtype host_arr[N];
  sunrealtype* some_arr = static_cast<sunrealtype*>(some_memory->ptr);
  if (mem_type == SUNMEMTYPE_DEVICE)
  {
    for (int i = 0; i < N; i++) { host_arr[i] = i * sunrealtype{1.0}; }
    myQueue.memcpy(some_memory->ptr, host_arr, bytes_to_alloc);
    myQueue.wait_and_throw();
    some_arr = host_arr;
  }
  else
  {
    for (int i = 0; i < N; i++) { some_arr[i] = i * sunrealtype{1.0}; }
  }
  if (print_test_status) std::cout << "  SUNMemoryHelper_Alloc... PASSED\n";

  // Try and copy the memory
  if (print_test_status) std::cout << "  SUNMemoryHelper_Copy... \n";
  SUNMemory other_memory = nullptr;
  retval = SUNMemoryHelper_Alloc(helper, &other_memory, bytes_to_alloc,
                                 mem_type, static_cast<void*>(&myQueue));
  if (retval)
  {
    if (print_test_status) std::cout << "  SUNMemoryHelper_Alloc... FAILED\n";
    return -1;
  }
  retval = SUNMemoryHelper_Copy(helper, other_memory, some_memory,
                                bytes_to_alloc, static_cast<void*>(&myQueue));
  if (retval)
  {
    if (print_test_status)
      std::cout << "  SUNMemoryHelper_Copy... FAILED retval\n";
    return -1;
  }
  if (mem_type == SUNMEMTYPE_DEVICE)
  {
    sunrealtype other_arr[N];
    myQueue.memcpy(other_arr, other_memory->ptr, bytes_to_alloc);
    myQueue.wait_and_throw();
    for (int i = 0; i < N; i++)
    {
      if (some_arr[i] != other_arr[i])
      {
        if (print_test_status)
          std::cout << "  SUNMemoryHelper_Copy... FAILED comparison\n";
        return -1;
      }
    }
  }
  else
  {
    sunrealtype* other_arr = static_cast<sunrealtype*>(other_memory->ptr);
    for (int i = 0; i < N; i++)
    {
      if (some_arr[i] != other_arr[i])
      {
        if (print_test_status)
          std::cout << "  SUNMemoryHelper_Copy... FAILED comparison\n";
        return -1;
      }
    }
  }
  if (print_test_status) std::cout << "  SUNMemoryHelper_Copy... PASSED\n";

  // Try and deallocate
  if (print_test_status) std::cout << "  SUNMemoryHelper_Dealloc... \n";
  retval = SUNMemoryHelper_Dealloc(helper, some_memory,
                                   static_cast<void*>(&myQueue));
  if (retval)
  {
    if (print_test_status) std::cout << "  SUNMemoryHelper_Dealloc... FAILED\n";
    return -1;
  }
  retval = SUNMemoryHelper_Dealloc(helper, other_memory,
                                   static_cast<void*>(&myQueue));
  if (retval)
  {
    if (print_test_status) std::cout << "  SUNMemoryHelper_Dealloc... FAILED\n";
    return -1;
  }
  if (print_test_status) std::cout << "  SUNMemoryHelper_Dealloc... PASSED\n";

  // Check alloc stats
  if (print_test_status) std::cout << "  SUNMemoryHelper_GetAllocStats... \n";
  unsigned long num_allocations, num_deallocations;
  size_t bytes_allocated, bytes_high_watermark;

  retval = SUNMemoryHelper_GetAllocStats(helper, mem_type, &num_allocations,
                                         &num_deallocations, &bytes_allocated,
                                         &bytes_high_watermark);
  if (retval)
  {
    if (print_test_status)
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    return -1;
  }
  if (print_test_status)
    std::cout << "\tnum_allocations = " << num_allocations
              << " num_deallocations = " << num_deallocations
              << " bytes_allocated = " << bytes_allocated
              << " bytes_high_watermark = " << bytes_high_watermark << "\n";
  if (num_allocations != 2)
  {
    if (print_test_status)
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    if (print_test_status) std::cout << "    num_allocations != 2\n";
    return -1;
  }
  if (num_deallocations != 2)
  {
    if (print_test_status)
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    if (print_test_status) std::cout << "    num_deallocations != 2\n";
    return -1;
  }
  if (bytes_allocated != 0)
  {
    if (print_test_status)
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    if (print_test_status) std::cout << "    bytes_allocated != 0\n";
    return -1;
  }
  if (bytes_high_watermark != bytes_to_alloc * 2)
  {
    if (print_test_status)
      std::cout << "  SUNMemoryHelper_GetAllocStats... FAILED\n";
    if (print_test_status) std::cout << "    bytes_high_watermark != 0\n";
    return -1;
  }
  if (print_test_status)
    std::cout << "  SUNMemoryHelper_GetAllocStats... PASSED\n";
  return retval;
}

int main(int argc, char* argv[])
{
  sundials::Context sunctx;

  std::cout << "Testing the SUNMemoryHelper_Sycl module... \n";

  std::cout << "  SUNMemoryHelper_Sycl... \n";
  SUNMemoryHelper helper = SUNMemoryHelper_Sycl(sunctx);
  if (!helper)
  {
    std::cout << "  SUNMemoryHelper_Sycl... FAILED\n";
    return -1;
  }
  std::cout << "  SUNMemoryHelper_Sycl... PASSED\n";

  std::cout << "With host memory... \n";
  test_instance(helper, SUNMEMTYPE_HOST, true);
  std::cout << "With pinned memory... \n";
  test_instance(helper, SUNMEMTYPE_PINNED, true);
  std::cout << "With device memory... \n";
  test_instance(helper, SUNMEMTYPE_DEVICE, true);
  std::cout << "With uvm memory... \n";
  test_instance(helper, SUNMEMTYPE_UVM, true);

  std::cout << "  SUNMemoryHelper_Clone... \n";
  SUNMemoryHelper helper2 = SUNMemoryHelper_Clone(helper);
  if (!helper || test_instance(helper2, SUNMEMTYPE_HOST, false))
  {
    std::cout << "  SUNMemoryHelper_Clone... FAILED\n";
    return -1;
  }
  std::cout << "  SUNMemoryHelper_Clone... PASSED\n";

  // Check destroy
  std::cout << "  SUNMemoryHelper_Destroy... \n";
  if (SUNMemoryHelper_Destroy(helper) || SUNMemoryHelper_Destroy(helper2))
  {
    std::cout << "  SUNMemoryHelper_Destroy... FAILED\n";
    return -1;
  }
  std::cout << "  SUNMemoryHelper_Destroy... PASSED\n";

  return 0;
}
