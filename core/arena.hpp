#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace inferno::core {

// Bump allocator over fixed-size blocks. Individual objects are never freed;
// the whole arena is released at once. Used for factorization scratch space
// and other short-lived, high-churn allocations where malloc/free overhead
// would dominate.
class Arena {
 public:
  explicit Arena(std::size_t block_size = 1u << 20) : block_size_(block_size) {}

  void* Allocate(std::size_t bytes, std::size_t alignment = alignof(std::max_align_t)) {
    if (blocks_.empty() || !FitsInCurrentBlock(bytes, alignment)) {
      AddBlock(std::max(bytes + alignment, block_size_));
    }
    Block& b = blocks_.back();
    std::size_t aligned = AlignUp(b.used, alignment);
    void* ptr = b.data.get() + aligned;
    b.used = aligned + bytes;
    return ptr;
  }

  template <typename T>
  T* AllocateArray(std::size_t count) {
    return static_cast<T*>(Allocate(sizeof(T) * count, alignof(T)));
  }

  void Reset() {
    for (auto& b : blocks_) b.used = 0;
  }

  std::size_t BlockCount() const { return blocks_.size(); }

 private:
  struct Block {
    std::unique_ptr<std::byte[]> data;
    std::size_t capacity = 0;
    std::size_t used = 0;
  };

  static std::size_t AlignUp(std::size_t offset, std::size_t alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
  }

  bool FitsInCurrentBlock(std::size_t bytes, std::size_t alignment) const {
    const Block& b = blocks_.back();
    std::size_t aligned = AlignUp(b.used, alignment);
    return aligned + bytes <= b.capacity;
  }

  void AddBlock(std::size_t min_size) {
    Block b;
    b.capacity = min_size;
    b.data = std::make_unique<std::byte[]>(min_size);
    b.used = 0;
    blocks_.push_back(std::move(b));
  }

  std::size_t block_size_;
  std::vector<Block> blocks_;
};

}  // namespace inferno::core
