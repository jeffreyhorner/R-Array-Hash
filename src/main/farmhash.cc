#include <farmhash.h>
extern "C" {
  size_t StringHash(const char* s, size_t len) {
    return util::Hash(s, len);
  }
}

