#include <iostream>
#include "floating.cpp"

#define ASSERT(expr, ...)   if(!(expr)){ \
                              printf(__VA_ARGS__); \
                            }
#define ASSERT_EQ(left, right, ...) \
  if((left) != (right)){ \
    std::cout << \
    __VA_ARGS__ << "\n" \
    "correct would be: " << (left) << "\n" \
    "but the result is:" << (right) << '\n' << std::endl; \
    std::exit(1); \
  }
#define ASSERT_ALMOST_EQ(left, right, ...) \
  if(!EQ(left,right)){ \
    std::cout << \
    __VA_ARGS__ << "\n" \
    "correct would be: " << (left) << "\n" \
    "but the result is:" << (right) << '\n' << std::endl; \
    std::exit(1); \
  }

int main() {
  // static_assert(!is_nan(Float32::INFTY));
  // static_assert(Float32::INFTY.bits!=Float32::NEG_INFTY.bits);
  constexpr float TEST_FLOATS[] = {
      0.1f,
      0.5f,
      1.75f,
      10.10f,
      0.001f,
      321325432.321f,
      -543254325.5432f,
      0.0000000413f,
      1e-20,
      std::numeric_limits<float>::infinity(),
      -std::numeric_limits<
          float>::infinity() /*, std::numeric_limits<float>::max() */,
      std::numeric_limits<float>::quiet_NaN(),
      /* std::numeric_limits<float>::min(), */ -0.0f,
      0.0f};

  // test pack/unpack
  for (float f : TEST_FLOATS) {
    Float32 my_float{std::bit_cast<uint32_t>(f)};
    Float32 repacked{pack(unpack(my_float))};

    ASSERT_EQ(my_float, repacked, "pack/unpack identity");
  }

  // test unary
  for (const float f : TEST_FLOATS) {
    const Float32 fa{std::bit_cast<uint32_t>(f)};

    // std::cout << "current unary argument: " << f << std::endl;

    ASSERT_EQ(fa, pack(unpack(fa)), "pack/unpack identity");

    ASSERT_EQ(Float32{std::bit_cast<uint32_t>(-f)}, neg(fa),
              "negation falied for: " << f);
    ASSERT_EQ(std::abs(f), abs(fa), "abs falied for: " << f);
    ASSERT_EQ(std::abs(f), abs(fa), "abs falied for: " << f);

    ASSERT_EQ(f == 0.0f, is_zero(fa), "is_zero falied for: " << f);
    ASSERT_EQ(std::isnan(f), is_nan(fa), "is_nan falied for: " << f);
    ASSERT_EQ(std::isinf(f), is_inf(fa), "is_inf falied for: " << f);
    ASSERT_EQ(std::signbit(f), is_sign_minus(fa), "signbit falied for: " << f);

    ASSERT_ALMOST_EQ(Float32{std::bit_cast<uint32_t>(std::sqrt(f))}, sqrt(fa),
                     "Wrong result of sqrt(" << f << ")");

    ASSERT_EQ(static_cast<int32_t>(f), round_to_int<int32_t>(fa),
              "to int32 conversion failed for: " << f);
    ASSERT_EQ(static_cast<float>(static_cast<int32_t>(f)),
              convert_from_int<int32_t>(round_to_int<int32_t>(fa)),
              "from int32 conversion failed for: " << static_cast<int32_t>(f));
  }

  /*
   test binary operations
  */
  for (auto a : TEST_FLOATS) {
    for (auto b : TEST_FLOATS) {
      Float32 af{std::bit_cast<uint32_t>(a)};
      Float32 bf{std::bit_cast<uint32_t>(b)};

      // COMPARISON
      ASSERT_EQ(a < b, qls(af, bf),
                "less than comparison failed for: " << a << "<" << b << '\n');
      ASSERT_EQ(a <= b, qle(af, bf),
                "less than or eq comparison failed for: " << a << "<=" << b
                                                          << '\n');
      ASSERT_EQ(a > b, qgr(af, bf),
                "greater than comparison failed for: " << a << ">" << b
                                                       << '\n');
      ASSERT_EQ(a >= b, qge(af, bf),
                "greater than or eq comparison failed for: " << a << ">=" << b
                                                             << '\n');
      ASSERT_EQ(a == b, qeq(af, bf),
                "eq comparison failed for: " << a << "==" << b << '\n');

      ASSERT_ALMOST_EQ(Float32{std::bit_cast<uint32_t>(a * b)}, mul(af, bf),
                       "Wrong result multiplication of " << a << " * " << b);
      ASSERT_ALMOST_EQ(Float32{std::bit_cast<uint32_t>(a + b)}, add(af, bf),
                       "Wrong result addition of " << a << " + " << b);
      ASSERT_ALMOST_EQ(Float32{std::bit_cast<uint32_t>(a - b)}, sub(af, bf),
                       "Wrong result subtraction of " << a << " - " << b);
      ASSERT_ALMOST_EQ(Float32{std::bit_cast<uint32_t>(a / b)}, div(af, bf),
                       "Wrong result division of " << a << " / " << b);
    }
  }
}
