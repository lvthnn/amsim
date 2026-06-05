#pragma once

#include <cstdint>

namespace amsim {

enum class Generation { Current, Parents };

enum class Sex : uint8_t { Unknown = 0, Male = 1, Female = 2 };

}  // namespace amsim
