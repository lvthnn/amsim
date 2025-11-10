#include <fstream>

#include <amsim/logger.h>
#include <amsim/log_level.h>

int main() {
  LOG_FILE("test.log", LogLevel::DEBUG);
  LOG_DEBUG("debug message");
  LOG_INFO("info message");
  LOG_WARNING("warning message");
  LOG_ERROR("error message");
}
