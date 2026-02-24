#pragma once

#include <amsim/init.h>
#include <h5pp/h5pp.h>

#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <unordered_set>

namespace amsim {

using Job = std::function<void()>;

class Writer {
 public:
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  static Writer& get_instance(const Simulation& sim) {
    static Writer instance(sim);
    return instance;
  }

  ~Writer() {
    {
      std::lock_guard<std::mutex> lg(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    if (thread_.joinable()) thread_.join();
  }

  template <typename T>
  void create(
      T&& data,
      const std::string& name,
      std::size_t rep,
      std::size_t n_rows,
      std::size_t n_cols,
      const std::vector<std::string>& row_labels,
      const std::vector<std::string>& col_labels);

  template <typename T>
  void write(
      T&& data, const std::string& name, std::size_t rep, std::size_t gen);

 private:
  explicit Writer(const Simulation& sim)
      : file_(sim.output_dir, h5pp::FileAccess::COLLISION_FAIL),
        n_gens_(sim.n_generations) {
    thread_ = std::thread(&Writer::threadCallback, this);
  }

  h5pp::File file_;
  std::size_t n_gens_;

  std::unordered_set<std::string> datasets_;
  std::deque<Job> jobs_;
  std::condition_variable cv_;
  std::thread thread_;
  std::mutex mutex_;
  bool done_ = false;

  void threadCallback();
};

inline void Writer::threadCallback() {
  std::unique_lock<std::mutex> lk(mutex_);
  while (!done_ || !jobs_.empty()) {
    cv_.wait(lk, [this] { return done_ || !jobs_.empty(); });

    while (!jobs_.empty()) {
      Job job = std::move(jobs_.front());
      jobs_.pop_front();
      lk.unlock();
      job();
      lk.lock();
    }
  }
}

template <typename T>
inline void Writer::create(
    T&& /*data*/,
    const std::string& name,
    std::size_t rep,
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels) {
  {
    std::lock_guard<std::mutex> lg(mutex_);

    if (!file_.attributeExists(name, "row_labels"))
      file_.writeAttribute(row_labels, name, "row_labels");

    if (!file_.attributeExists(name, "col_labels"))
      file_.writeAttribute(col_labels, name, "col_labels");

    std::string name_rep = std::format("{}/rep{:03d}", name, rep);

    if (!datasets_.contains(name_rep)) {
      file_.createDataset<T>(
          T{},
          name_rep,
          std::vector<hsize_t>{n_gens_, n_rows, n_cols},
          H5D_CHUNKED);

      datasets_.insert(name_rep);
    }
  }
}

template <typename T>
inline void Writer::write(
    T&& data, const std::string& name, std::size_t rep, std::size_t gen) {
  {
    std::lock_guard<std::mutex> lg(mutex_);

    jobs_.push_back([this, copy = data, path])
  }
}


}  // namespace amsim
