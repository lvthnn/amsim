// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <boost/asio.hpp>
#include <boost/process.hpp>
#include <format>
#include <semaphore>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace amsim::process {

namespace bp = boost::process;
namespace asio = boost::asio;

struct ProcessResult {
  int exit_code;
  std::string output;
};

inline std::counting_semaphore<64> process_bound{static_cast<ptrdiff_t>(
    std::thread::hardware_concurrency() > 0
        ? std::thread::hardware_concurrency()
        : 4)};

// resolve an executable to be run with bp::process
inline bp::filesystem::path resolveExecutable(const std::string& s) {
  bp::filesystem::path exec(s);

  // if the exec is similar to a filesystem path, the user is directing us to
  // a specific path — in which case, resolve it and make sure it exists
  if (exec.has_parent_path()) {
    bp::filesystem::path exec_abs = bp::filesystem::absolute(exec);
    if (!bp::filesystem::exists(exec_abs))
      throw std::runtime_error(
          std::format("Executable {} not found", exec_abs.string()));

    return exec_abs;
  }

  // otherwise, we should search PATH for the executable
  auto found = bp::environment::find_executable(exec);
  if (found.empty())
    throw std::runtime_error(
        std::format("Executable {} not found in PATH", exec.string()));

  return found;
}

// very similar to resolveExecutable, but just a premature check — is it
// resolvable?
inline void checkProcessAvailable(const std::string& s) {
  bp::filesystem::path exec(s);

  if (exec.has_parent_path()) {
    bp::filesystem::path exec_abs = bp::filesystem::absolute(exec_abs);
    if (bp::filesystem::exists(exec_abs)) {
      throw std::runtime_error(
          std::format("Executable {} not found", exec_abs.string()));
    }
  }

  auto found = bp::environment::find_executable(exec);
  if (found.empty())
    throw std::runtime_error(
        std::format("Executable {} not found in PATH", exec.string()));
}

// Spawns a process and runs a command with arguments, capturing output from
// stderr and stdout, and returns a ProcessResult
inline ProcessResult runProcess(
    const std::string& cmd, const std::vector<std::string>& args) {
  // declare I/O context, a pipe to read stdout + stderr output from,
  // and a string to capture the output of the two
  asio::io_context ctx;
  asio::readable_pipe pipe{ctx};
  std::string pipe_out;

  // resolve the executable
  bp::filesystem::path exec = resolveExecutable(cmd);

  // launch the process
  bp::process proc(
      ctx, exec, args, bp::process_stdio{.out = pipe, .err = pipe});

  // read in stdout + stderr
  boost::system::error_code ec;
  asio::read(pipe, asio::dynamic_buffer(pipe_out), ec);

  // the error code is non-trivial — is it something different from the excepted
  // asio::error::eof?
  if (ec && (ec != asio::error::eof))
    throw std::runtime_error(
        std::format(
            "process {} failed to parse stdout / stderr output",
            exec.string()));

  proc.wait();

  return ProcessResult{.exit_code = proc.exit_code(), .output = pipe_out};
}

inline ProcessResult runProcessBounded(
    const std::string& cmd, const std::vector<std::string>& args) {
  process_bound.acquire();
  ProcessResult result = runProcess(cmd, args);
  process_bound.release();
  return result;
}

}  // namespace amsim::process
