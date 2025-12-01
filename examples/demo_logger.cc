// This example demonstrates the functionality of the simulation logger
#include <amsim/log_level.h>
#include <amsim/logger.h>

#include <fstream>

int main() {
  LOG_FILE("test.log", amsim::LogLevel::DEBUG);
  LOG_DEBUG("debug message");
  LOG_INFO("info message");
  LOG_WARNING("warning message");
  LOG_ERROR("error message");
}
