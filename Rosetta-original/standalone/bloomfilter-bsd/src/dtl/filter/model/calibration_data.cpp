#include "calibration_data.hpp"

#include <sstream>

namespace dtl {
namespace filter {
namespace model {

using bbf_config_t = dtl::blocked_bloomfilter_config;
using cf_config_t = dtl::cuckoofilter::config;
using timings_t = std::vector<timing>;

//===----------------------------------------------------------------------===//
/// Open and map the calibration data file.
void
calibration_data::open_file() {
  assert(file_descriptor_ == -1);
  assert(file_size_ == 0);
  assert(mapped_data_ == nullptr);

  // test if the file exists
  file_descriptor_ = open(filename_.c_str(), O_RDONLY);
  if (file_descriptor_ != -1) {
    // file exists and is ready to be read

    // drop in-memory state
    cache_sizes_.clear();
    filter_sizes_.clear();
    bbf_tuning_params_.clear();
    cf_tuning_params_.clear();
    bbf_delta_tls_.clear();
    cf_delta_tls_.clear();
    skyline_matrix_mem_.clear();
    changed_ = false;

    struct stat sb;
    fstat(file_descriptor_, &sb);
    file_size_ = static_cast<u64>(sb.st_size);

    if (file_size_ < 24) {
      throw std::runtime_error("Failed to read file '" + filename_ + "'. Illegal format.");
    }

    mapped_data_ = reinterpret_cast<$u8*>(mmap(nullptr, file_size_, PROT_READ, MAP_PRIVATE, file_descriptor_, 0));
    if (mapped_data_ == nullptr) {
      close(file_descriptor_);
      throw std::runtime_error("Failed to mmap file '" + filename_ + "'.");
    }

    // check signature
    if (std::memcmp(file_signature, mapped_data_, file_signature_size) != 0) {
      munmap(mapped_data_, file_size_);
      close(file_descriptor_);
      throw std::runtime_error("Failed to read file '" + filename_ + "'. Illegal file type.");
    }

    // check version
    if (reinterpret_cast<u64*>(mapped_data_)[1] != file_version) {
      munmap(mapped_data_, file_size_);
      close(file_descriptor_);
      throw std::runtime_error("Failed to read file '" + filename_ + "'. Version mismatch.");
    }

    // read the header
    $u64* header = reinterpret_cast<$u64*>(mapped_data_ + 16);
    auto mem_levels = *header;
    header++;
    cache_sizes_.clear();
    for ($u64 i = 0; i < mem_levels - 1; i++) {
      cache_sizes_.push_back(*header);
      header++;
    }
    filter_sizes_.clear();
    for ($u64 i = 0; i < mem_levels; i++) {
      filter_sizes_.push_back(*header);
      header++;
    }

    // read the number of bloom filter configurations
    u64 bbf_config_count = *header;
    header++;
    // read the number of cuckoo filter configurations
    u64 cf_config_count = *header;

    // calculate the header size
    std::size_t header_size = file_signature_size
        + sizeof(u64) // version
        + sizeof(u64) // mem levels
        + (mem_levels - 1) * sizeof(u64) // cache sizes
        + mem_levels * sizeof(u64) // filter sizes
        + sizeof(u64) // number of bloom filter configurations
        + sizeof(u64) // number of cuckoo filter configurations
    ;

    // init begin and end pointers
    bbf_config_begin_ = reinterpret_cast<bbf_config_t*>(&mapped_data_[header_size]);
    bbf_config_end_ = bbf_config_begin_ + bbf_config_count;
    bbf_timing_begin_ = reinterpret_cast<timing*>(bbf_config_end_);
    bbf_timing_end_ = bbf_timing_begin_ + mem_levels * bbf_config_count;
    bbf_tuning_params_begin_ = reinterpret_cast<tuning_params*>(bbf_timing_end_);
    bbf_tuning_params_end_ = bbf_tuning_params_begin_ + bbf_config_count;
    cf_config_begin_ = reinterpret_cast<cf_config_t*>(bbf_tuning_params_end_);
    cf_config_end_ = cf_config_begin_ + cf_config_count;
    cf_timing_begin_ = reinterpret_cast<timing*>(cf_config_end_);
    cf_timing_end_ = cf_timing_begin_ + mem_levels * cf_config_count;
    cf_tuning_params_begin_ = reinterpret_cast<tuning_params*>(cf_timing_end_);
    cf_tuning_params_end_ = cf_tuning_params_begin_ + cf_config_count;

    // check if the file contains a skyline matrix
    const auto file_pos = std::distance(mapped_data_, reinterpret_cast<$u8*>(cf_tuning_params_end_));
    if (file_pos < file_size_) {
      skyline_matrix_ = reinterpret_cast<skyline_matrix*>(cf_tuning_params_end_);
    }
    else {
      init_skyline_matrix();
    }

    changed_ = false;
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Unmap and close the calibration data file.
void
calibration_data::close_file() {
  if (mapped_data_) {
    munmap(mapped_data_, file_size_);
  }
  if (file_descriptor_ != -1) {
    close(file_descriptor_);
  }
  mapped_data_ = nullptr;
  file_descriptor_ = -1;
  file_size_ = 0;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Write the changes to the file.
void
calibration_data::persist() {
  if (! changed()) return;
  auto serialized_data = serialize();
  // unmap the existing file
  close_file();

  // write a new file
  std::ofstream file(filename_, std::ios::binary | std::ios::out);
  file.write((char*) &serialized_data[0], serialized_data.size());
  file.close();

  // map the newly created file
  open_file();
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Serialize the calibration data.
std::vector<$u8>
calibration_data::serialize() {

  // read data from file to memory
  u64 mem_levels = get_mem_levels();
  {
    // bloom
    u64 cnt = std::distance(bbf_config_begin_, bbf_config_end_);
    bbf_config_t* config_reader = bbf_config_begin_;
    timing* timing_reader = bbf_timing_begin_;
    tuning_params* tuning_reader = bbf_tuning_params_begin_;
    for (std::size_t i = 0; i < cnt; i++) {
      // check if already in memory
      auto search_delta_tls = bbf_delta_tls_.find(*config_reader);
      if (search_delta_tls == bbf_delta_tls_.end()) {
        // read from file
        timings_t delta_timings;
        for (std::size_t j = 0; j < mem_levels; j++) {
          delta_timings.push_back(*timing_reader);
          timing_reader++;
        }
        bbf_delta_tls_.insert(std::make_pair(*config_reader, delta_timings));
      }
      else {
        timing_reader += mem_levels;
      }
      auto search_tuning_params = bbf_tuning_params_.find(*config_reader);
      if (search_tuning_params == bbf_tuning_params_.end()) {
        bbf_tuning_params_.insert(std::make_pair(*config_reader, *tuning_reader));
      }
      config_reader++;
      tuning_reader++;
    }
  }
  {
    // cuckoo
    u64 cnt = std::distance(cf_config_begin_, cf_config_end_);
    cf_config_t* config_reader = cf_config_begin_;
    timing* timing_reader = cf_timing_begin_;
    tuning_params* tuning_reader = bbf_tuning_params_begin_;
    for (std::size_t i = 0; i < cnt; i++) {
      // check if already in memory
      auto search_timings = cf_delta_tls_.find(*config_reader);
      if (search_timings == cf_delta_tls_.end()) {
        // read from file
        timings_t delta_timings;
        for (std::size_t j = 0; j < mem_levels; j++) {
          delta_timings.push_back(*timing_reader);
          timing_reader++;
        }
        cf_delta_tls_.insert(std::make_pair(*config_reader, delta_timings));
      }
      else {
        timing_reader += mem_levels;
      }
      auto search_tuning_params = cf_tuning_params_.find(*config_reader);
      if (search_tuning_params == cf_tuning_params_.end()) {
        cf_tuning_params_.insert(std::make_pair(*config_reader, *tuning_reader));
      }
      config_reader++;
      tuning_reader++;
    }
  }

  u64 bbf_config_count = bbf_delta_tls_.size();
  u64 cf_config_count = cf_delta_tls_.size();

  std::vector<$u8> buffer;

  // calculate the header size
  const std::size_t header_size = file_signature_size
      + sizeof(u64) // version
      + sizeof(u64) // mem levels
      + (get_mem_levels() - 1) * sizeof(u64) // cache sizes
      + get_mem_levels() * sizeof(u64) // filter sizes
      + sizeof(u64) // number of bloom filter configurations
      + sizeof(u64) // number of cuckoo filter configurations
  ;

  const std::size_t buffer_size = header_size
      + bbf_config_count * (sizeof(bbf_config_t) + get_mem_levels() * sizeof(timing) + sizeof(tuning_params)) // bbf config + timings + tuning parameters
      + cf_config_count * (sizeof(cf_config_t) + get_mem_levels() * sizeof(timing) + sizeof(tuning_params)) // cf config + timings + tuning parameters
      + get_skyline_matrix_ptr()->size_in_bytes() // skyline matrix
  ;

  // allocate a buffer
  buffer.resize(buffer_size);

  {
    // write header
    auto* writer = &buffer[0];

    // file signature
    std::memcpy(writer, file_signature, file_signature_size);
    writer += file_signature_size;

    *reinterpret_cast<$u64*>(writer) = file_version;
    writer += sizeof(u64);

    // mem levels
    *reinterpret_cast<$u64*>(writer) = get_mem_levels();
    writer += sizeof(u64);

    // cache sizes
    for (u64 cache_size : cache_sizes_) {
      *reinterpret_cast<$u64*>(writer) = cache_size;
      writer += sizeof(u64);
    }

    // filter sizes
    for (u64 filter_size : filter_sizes_) {
      *reinterpret_cast<$u64*>(writer) = filter_size;
      writer += sizeof(u64);
    }
  }

  {

    // number of configurations
    *reinterpret_cast<$u64*>(&buffer[0] + header_size - 2 * sizeof(u64)) = bbf_config_count;
    *reinterpret_cast<$u64*>(&buffer[0] + header_size - 1 * sizeof(u64)) = cf_config_count;

    // configurations, timings, and tuning params
    u64 bbf_config_offset = header_size;
    u64 bbf_timing_offset = bbf_config_offset + bbf_config_count * sizeof(bbf_config_t);
    u64 bbf_tuning_offset = bbf_timing_offset + bbf_config_count * get_mem_levels() * sizeof(timing);
    bbf_config_t* bf_configs = reinterpret_cast<bbf_config_t*>(&buffer[0] + bbf_config_offset);
    timing* bf_timings = reinterpret_cast<timing*>(&buffer[0] + bbf_timing_offset);
    tuning_params* bf_tuning = reinterpret_cast<tuning_params*>(&buffer[0] + bbf_tuning_offset);
    for (auto& ct : bbf_delta_tls_) {
      *bf_configs = ct.first;
      bf_configs++;
      for (auto& t : ct.second) {
        *bf_timings = t;
        bf_timings++;
      }
      auto search = bbf_tuning_params_.find(ct.first);
      *bf_tuning = search->second;
      bf_tuning++;
    }
    u64 cf_config_offset = bbf_tuning_offset + bbf_config_count * sizeof(tuning_params);
    u64 cf_timing_offset = cf_config_offset + cf_config_count * sizeof(cf_config_t);
    u64 cf_tuning_offset = cf_timing_offset + cf_config_count * get_mem_levels() * sizeof(timing);
    cf_config_t* cf_configs = reinterpret_cast<cf_config_t*>(&buffer[0] + cf_config_offset);
    timing* cf_timings = reinterpret_cast<timing*>(&buffer[0] + cf_timing_offset);
    tuning_params* cf_tuning = reinterpret_cast<tuning_params*>(&buffer[0] + cf_tuning_offset);
    for (auto& ct : cf_delta_tls_) {
      *cf_configs = ct.first;
      cf_configs++;
      for (auto& t : ct.second) {
        *cf_timings = t;
        cf_timings++;
      }
      auto search = cf_tuning_params_.find(ct.first);
      *cf_tuning = search->second;
      cf_tuning++;
    }

    u64 bf_skyline_offset = cf_tuning_offset + cf_config_count * sizeof(tuning_params);
    get_skyline_matrix_ptr()->serialize(&buffer[bf_skyline_offset]);

  }

  return buffer;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Add the delta-t_l values for the given filter configuration.
/// Note: All changes are transient, until the persist function is called.
void
calibration_data::put_timings(const bbf_config_t& config, const timings_t& delta_timings) {
  // ensure, that we have tuning parameters for the given config
  tuning_params params;
  $u1 add_tuning_params = false;
  try {
    params = get_tuning_params(config);
  } catch (...) {
    // no timings found for the given configuration
    params = get_null_tuning_params();
    add_tuning_params = true;
  }

  bbf_delta_tls_.insert(std::make_pair(config, delta_timings));

  if (add_tuning_params) {
    put_tuning_params(config, params);
  }
  changed_ = true;
}

void
calibration_data::put_timings(const cf_config_t& config, const timings_t& delta_timings) {
  // ensure, that we have tuning parameters for the given config
  tuning_params params;
  $u1 add_tuning_params = false;
  try {
    params = get_tuning_params(config);
  } catch (...) {
    // no timings found for the given configuration
    params = get_null_tuning_params();
    add_tuning_params = true;
  }

  cf_delta_tls_.insert(std::make_pair(config, delta_timings));

  if (add_tuning_params) {
    put_tuning_params(config, params);
  }
  changed_ = true;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Get the delta-t_l values for the given filter configuration.
timings_t
calibration_data::get_timings(const bbf_config_t& config) const {
  auto search = bbf_delta_tls_.find(config);
  if (search != bbf_delta_tls_.end()) {
    return search->second;
  }
  // search in the file
  const auto found = std::lower_bound(bbf_config_begin_, bbf_config_end_, config);
  if (!found || config != *found) {
    throw std::runtime_error("Failed to find configuration.");
  }
  const auto config_idx = std::distance(bbf_config_begin_, found);

  // read the timings
  timings_t timings;
  for (std::size_t i = 0; i < get_mem_levels(); i++) {
    timings.push_back(bbf_timing_begin_[get_mem_levels() * config_idx + i]);
  }
  return timings;
}

u1
calibration_data::has_timings(const bbf_config_t& config) const {
  const auto delta_timings = get_timings(config);
  if (delta_timings.empty() || delta_timings[0].cycles_per_lookup == 0.0 || delta_timings[0].nanos_per_lookup == 0.0) {
    return false;
  }
  return true;
}


timings_t
calibration_data::get_timings(const cf_config_t& config) const {
  auto search = cf_delta_tls_.find(config);
  if (search != cf_delta_tls_.end()) {
    return search->second;
  }
  // search in the file
  const auto found = std::lower_bound(cf_config_begin_, cf_config_end_, config);
  if (!found || config != *found) {
    throw std::runtime_error("Failed to find configuration.");
  }
  const auto config_idx = std::distance(cf_config_begin_, found);

  // read the timings
  timings_t timings;
  for (std::size_t i = 0; i < get_mem_levels(); i++) {
    timings.push_back(cf_timing_begin_[get_mem_levels() * config_idx + i]);
  }
  return timings;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Add the tuning parameters for the given filter configuration.
/// Note: All changes are transient, until the persist function is called.
void
calibration_data::put_tuning_params(const bbf_config_t& config, const tuning_params& params) {
  // ensure, that we have timings for the given config
  timings_t timings;
  $u1 add_timings = false;
  try {
    timings = get_timings(config);
  } catch (...) {
    // no timings found for the given configuration
    timings = get_null_timings();
    add_timings = true;
  }

  bbf_tuning_params_.insert(std::make_pair(config, params));

  if (add_timings) {
    put_timings(config, timings);
  }
  changed_ = true;
}

void
calibration_data::put_tuning_params(const cf_config_t& config, const tuning_params& params) {
  // ensure, that we have timings for the given config
  timings_t timings;
  $u1 add_timings = false;
  try {
    timings = get_timings(config);
  } catch (...) {
    // no timings found for the given configuration
    timings = get_null_timings();
    add_timings = true;
  }

  cf_tuning_params_.insert(std::make_pair(config, params));

  if (add_timings) {
    put_timings(config, timings);
  }
  changed_ = true;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Get the tuning parameters for the given filter configuration.
tuning_params
calibration_data::get_tuning_params(const bbf_config_t& config) const {
  auto search = bbf_tuning_params_.find(config);
  if (search != bbf_tuning_params_.end()) {
    return search->second;
  }
  // search in the file
  const auto found = std::lower_bound(bbf_config_begin_, bbf_config_end_, config);
  if (!found || config != *found) {
    throw std::runtime_error("Failed to find configuration.");
  }
  const auto config_idx = std::distance(bbf_config_begin_, found);

  // read the timings
  tuning_params params = bbf_tuning_params_begin_[config_idx];
  return params;
}
tuning_params
calibration_data::get_tuning_params(const cf_config_t& config) const {
  auto search = cf_tuning_params_.find(config);
  if (search != cf_tuning_params_.end()) {
    return search->second;
  }
  // search in the file
  const auto found = std::lower_bound(cf_config_begin_, cf_config_end_, config);
  if (!found || config != *found) {
    throw std::runtime_error("Failed to find configuration.");
  }
  const auto config_idx = std::distance(cf_config_begin_, found);

  // read the timings
  tuning_params params = cf_tuning_params_begin_[config_idx];
  return params;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Returns all BBF configuration for which calibration data is available.
std::set<bbf_config_t>
calibration_data::get_bbf_configs() const {
  std::set<bbf_config_t> ret_val;
  for (auto* it = bbf_config_begin_; it != bbf_config_end_; it++) {
    if (has_timings(*it)) {
      ret_val.insert(*it);
    }
  }
  return ret_val;
}

/// Returns all CF configuration for which calibration data is available.
std::set<cf_config_t>
calibration_data::get_cf_configs() const {
  std::set<cf_config_t> ret_val;
  for (auto* it = cf_config_begin_; it != cf_config_end_; it++) {
    ret_val.insert(*it);
  }
  for (auto& it : cf_delta_tls_) {
    ret_val.insert(it.first);
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
void
calibration_data::print(std::ostream& os) const {
  std::stringstream str;
  const auto mem_levels = get_mem_levels();
  str << "Cache sizes:" << std::endl;
  const auto cache_sizes = get_cache_sizes();
  for (std::size_t i = 0; i < cache_sizes.size(); i++) {
    str << "L" << (i + 1) << ": " << cache_sizes[i] << std::endl;
  }
  const auto filter_sizes = get_filter_sizes();
  str << "Filter sizes:" << std::endl;
  for (std::size_t i = 0; i < filter_sizes.size(); i++) {
    str << "Level " << (i + 1) << ": " << filter_sizes[i] << std::endl;
  }

  for (auto* it = bbf_config_begin_; it != bbf_config_end_; it++) {
    const auto idx = std::distance(bbf_config_begin_, it);
    str << "config=[" << *it << "]"
        << ", tuning_params=[" << bbf_tuning_params_begin_[idx] << "]"
        << ", delta_timings=[" << bbf_timing_begin_[mem_levels * idx];
    for (std::size_t l = 1; l < mem_levels; l++) {
      str << "," << bbf_timing_begin_[mem_levels * idx + l];
    }
    str << "]" << std::endl;
  }

  auto& skyline_meta = get_skyline_matrix()->meta_data_;
  str << "Skyline matrix contains " << (skyline_meta.n_values_count_ * skyline_meta.tw_values_count_) << " entries." << std::endl;

  // TODO dump CF data

  os << str.str();
}
//===----------------------------------------------------------------------===//


} // namespace model
} // namespace filter
} // namespace dtl
