#pragma once

#include <condition_variable>
#include <deque>
#include <filesystem>
#include <format>
#include <functional>
#include <highfive/eigen.hpp>
#include <highfive/highfive.hpp>
#include <mutex>
#include <thread>
#include <unordered_set>

namespace amsim {

using Job = std::function<void()>;

class Writer {
 public:
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  static Writer& get_instance(
      const std::filesystem::path& out_dir, std::size_t n_gens) {
    static Writer instance(out_dir, n_gens);
    instance_ = &instance;
    return *instance_;
  }

  static Writer& get_instance() {
    if (instance_ == nullptr)
      throw std::runtime_error("Must initialise Writer before calling");
    return *instance_;
  }

  template <typename T>
  static void create(
      T&& data,
      const std::string& name,
      const std::vector<std::string>& row_labels,
      const std::vector<std::string>& col_labels);

  template <typename T>
  static void write(
      T&& data, const std::string& name, std::size_t rep, std::size_t gen);

  ~Writer() {
    {
      std::lock_guard<std::mutex> lg(mutex_);
      done_ = true;
    }
    cv_.notify_one();
    if (thread_.joinable()) thread_.join();
  }

 private:
  explicit Writer(const std::filesystem::path& out_dir, std::size_t n_gens)
      : file_(out_dir / "results.h5", HighFive::File::Overwrite),
        n_gens_(n_gens) {
    thread_ = std::thread(&Writer::threadCallback, this);
  }

  template <typename T>
  void createEstimator(
      T&& data,
      const std::string& name,
      const std::vector<std::string>& row_labels,
      const std::vector<std::string>& col_labels);

  template <typename T>
  void writeEstimator(
      T&& data, const std::string& name, std::size_t rep, std::size_t gen);

  inline static Writer* instance_ = nullptr;

  HighFive::File file_;
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
inline void Writer::createEstimator(
    T&& /*data*/,
    const std::string& name,
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels) {
  {
    std::lock_guard<std::mutex> lg(mutex_);

    jobs_.push_back(
        [this, name, row_copy = row_labels, col_copy = col_labels]() {
          if (!file_.exist(name)) file_.createGroup(name);

          HighFive::Group estimator = file_.getGroup(name);

          if (!estimator.hasAttribute("row_labels"))
            estimator.createAttribute("row_labels", row_copy);

          if (!estimator.hasAttribute("col_labels"))
            estimator.createAttribute("col_labels", col_copy);
        });
  }
  cv_.notify_one();
}

template <typename T>
inline void Writer::writeEstimator(
    T&& data, const std::string& name, std::size_t rep, std::size_t gen) {
  {
    std::lock_guard<std::mutex> lg(mutex_);

    jobs_.push_back(
        [this,
         copy = data,
         dir = std::format("{}/rep{:03d}", name, rep),
         path = std::format("{}/rep{:03d}/gen{:03d}", name, rep, gen)] {
          if (!file_.exist(dir)) file_.createGroup(dir);
          file_.createDataSet(path, copy);
        });
  }
  cv_.notify_one();
}

template <typename T>
inline void Writer::create(
    T&& data,
    const std::string& name,
    const std::vector<std::string>& row_labels,
    const std::vector<std::string>& col_labels) {
  Writer::get_instance().createEstimator(data, name, row_labels, col_labels);
}

template <typename T>
inline void Writer::write(
    T&& data, const std::string& name, std::size_t rep, std::size_t gen) {
  Writer::get_instance().writeEstimator(data, name, rep, gen);
}

}  // namespace amsim
