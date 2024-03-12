/*
  PRNG engine interface
 */

#ifndef _SRC_LIB_UTILS_PRNG_ENGINE_H
#define _SRC_LIB_UTILS_PRNG_ENGINE_H

#include <limits>

namespace lbcrypto {

/**
 * @brief Defines the inteface for PRNG engines used by OpenFHE. 
 * New PRNG engines can be added by implementing concrete classes 
 * inhereting from this interface.
 */
struct PrngEngine {
  
  using result_type = uint32_t;

  /**
   * @brief minimum value used by C+11 distribution generators when no lower
   * bound is explicitly specified by the user
   */
  static constexpr result_type min() {
    return std::numeric_limits<result_type>::min();
  }

  /**
   * @brief maximum value used by C+11 distribution generators when no upper
   * bound is explicitly specified by the user
   */
  static constexpr result_type max() {
    return std::numeric_limits<result_type>::max();
  }

  /**
   * @brief main call to the PRNG
   */
  virtual result_type operator()() = 0;

  virtual ~PrngEngine() = default;
};


struct PrngFactory {
  // all C++11 distributions used in OpenFHE work by default with uint32_t
  // a different data type can be specified if needed for a particular
  // architecture
  using result_type = PrngEngine::result_type;
  virtual std::shared_ptr<PrngEngine> make(result_type seed) = 0;
  virtual std::shared_ptr<PrngEngine> make(const std::array<result_type, 16>& seed) = 0;
  virtual std::shared_ptr<PrngEngine> make(const std::array<result_type, 16>& seed, result_type counter) = 0;
  virtual ~PrngFactory() = default;
};


}  // namespace lbcrypto

#endif
