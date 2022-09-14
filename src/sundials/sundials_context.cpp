#include "sundials/sundials_context.hpp"
#include "sundials_context_impl.h"

namespace sundials {

Context::Context(void* comm)
{
  sunctx = underlying_ptr_.get();
  SUNContext_Create(comm, &sunctx);
}

Context::~Context()
{
  SUNContext_Free(&sunctx);
}

} // namespace sundials
