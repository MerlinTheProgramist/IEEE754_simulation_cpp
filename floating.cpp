#include <atomic>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>

#include "floating.hpp"

#define ASSERT(expr, ...)   if(!(expr)){ \
                              printf(__VA_ARGS__); \
                            }
#define ASSERT_EQ(left, right, ...) if((left)!=(right)){ \
                                      std::cout << (left) << " != " << (right) << '\n'; \
                                      std::cout << __VA_ARGS__ << std::endl; \
                                    }

// enum class IEEE754Type{
//   SignalingNan,
//   QuietNan,
//   NegativeInfinity,
//   NegativeZero,
//   PositiveInfinity,
//   PositiveZero,
//   NegativeNormalized,
//   PositiveNormalized,
//   Denormalized
// };




const Float32 
    ZERO{Float32{0}},
    INFTY{pack(IEEEFields{false, bit_mask<int32_t>(Float32::EXPONENT_SIZE), 0})},
    NEG_INFTY{neg(INFTY)},
    NaN{pack(IEEEFields{false, bit_mask<int32_t>(Float32::EXPONENT_SIZE), 1})};


constexpr bool is_nan(const Float32 f){
  auto a = unpack(f);
  return a.exponent==bit_mask<int32_t>(Float32::EXPONENT_SIZE) && a.mantissa != 0;
}
// bool is_sNaN(Float32 f){
//   auto a = unpack(f);
//   return a.exponent==bit_mask<uint32_t>(EXPONENT_SIZE) && a.mantissa != 0 && ;
// }
constexpr bool is_inf(Float32 f){
 auto a = unpack(f);
  return a.exponent==bit_mask<int32_t>(Float32::EXPONENT_SIZE) && a.mantissa == 0;
}
constexpr bool is_zero(Float32 f){
  return unpack(f).mantissa==0;
}

constexpr Float32 int_to_float(int32_t n){
  if(n==0)
    return Float32::ZERO;

  uint32_t sign = 0;
  if(n<0){
    sign = 1<<Float32::SIGN_OFFSET;
    n = -n;
  }

  uint32_t mantissa = n;
  uint32_t exponent = Float32::EXPONENT_BIAS + Float32::MANTISSA_SIZE;

  while(mantissa > bit_mask<uint32_t>(sizeof(uint32_t)*8)){
    mantissa >>= 1;
    exponent += 1;
  }
  while(mantissa < bit_mask<uint32_t>(Float32::MANTISSA_SIZE)){
    mantissa <<= 1;
    exponent -= 1;
  }
  
  exponent <<= Float32::EXPONENT_OFFSET;
  mantissa &= bit_mask<uint32_t>(Float32::MANTISSA_SIZE);

  return Float32{sign | exponent | mantissa};
}

constexpr Float32 add(const Float32 fa, const Float32 fb){
    if(is_nan(fa))
      return fa;
    if(is_nan(fb))
      return fb;

    if(is_sign_minus(fa))
      return sub(fb, neg(fa));
    if(is_sign_minus(fb))
      return sub(fa, neg(fb));

    // they are both positive:
    if(is_inf(fa) || is_inf(fb))
      return Float32::MAX;
}


/*
  1.0011b * 101.001b
  1.1875 * 5.125 = 6.05
  1.001101 + 2^2*1.01001 = 
  {e:0 m:0011} * {e:2 m:01001} = {e:2 m:01100}
*/
  constexpr Float32 mul(Float32 fa,Float32 fb){
    if(is_nan(fa))
      return fa;
    if(is_nan(fb))
      return fb;

    if((is_inf(fa) && is_zero(fb)) 
    || (is_zero(fa) && is_inf(fb)))
      return NaN;

    // sign remains only if sign is diffrent
    const bool sign_bit = (fa.bits ^ fb.bits) >> Float32::SIGN_OFFSET;
    if(is_inf(fa) || is_inf(fb))
      return Float32{(static_cast<uint32_t>(sign_bit) << Float32::SIGN_OFFSET) | INFTY.bits};
    
    IEEEFields a = unpack(fa);
    IEEEFields b = unpack(fb);
        
    /*
      exponent are added
      mantissas are mulilied and lower bits are discarded via shift
    */
      
    int16_t exponent = a.exponent + b.exponent - Float32::BIAS;
    uint64_t mantissa = static_cast<uint64_t>(a.mantissa)
                      * static_cast<uint64_t>(b.mantissa);

    mantissa >>= Float32::MANTISSA_SIZE;

    // normalization
    if(mantissa != 0){
      // check if mantissa is greater than its maximum
      while(mantissa > bit_mask<uint64_t>(Float32::MANTISSA_SIZE)){
        exponent += 1;
        mantissa >>= 1;
      }
      // check if mantissa is less than its minimum
      while(mantissa < 0x800000){
        exponent -= 1;
        mantissa <<= 1;
      }
    }

    const Float32 result = pack(IEEEFields{sign_bit,exponent,static_cast<uint32_t>(mantissa)});

    // prevent nan
    if(is_nan(result))
      return sign_bit ? Float32::MIN : Float32::MAX;

    return result;
}



constexpr Float32 scaleB(const Float32 x, const int32_t N){
  if(N==0 || is_inf(x) || is_nan(x) || is_zero(x))
    return x;

  if(N >= Float32::BIAS)
    return Float32::MAX;
  if(N <= -Float32::BIAS)
    return Float32::MIN;
  
  auto a = unpack(x);
  a.exponent += N;
  auto f = pack(a);
  return is_inf(f) ? Float32::MAX : f;
}

constexpr Float32 neg(const Float32 x){
  if(is_nan(x))
    return x;
  return Float32{x.bits ^ (1<<(Float32::SIGN_OFFSET))};
}
constexpr Float32 abs(const Float32 x){
  if(is_nan(x))
    return x;
  return Float32{x.bits & bit_mask<uint32_t>(Float32::TOTAL_SIZE-1)};
}
constexpr Float32 copy_sign(const Float32 x, const Float32 y){
  if(is_nan(x))
    return x;
  return Float32{(x.bits & bit_mask<uint32_t>(Float32::SIGN_OFFSET)) | (y.bits & (1<<(Float32::TOTAL_SIZE-1))) };
}

constexpr bool qeq(const Float32 a , const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;
  if(is_zero(a) && is_zero(b))
    return true;
  return a.bits == b.bits;
}
constexpr bool qne(const Float32 a, const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;
  if(is_zero(a) && is_zero(b))
    return true;
  return a.bits == b.bits;
}
constexpr bool qgr(const Float32 a, const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;

  const bool sign_a = a.bits >> Float32::SIGN_OFFSET;
  const bool sign_b = b.bits >> Float32::SIGN_OFFSET;

  // both positive
  if(!sign_a && !sign_b)
    return a.bits > b.bits;
  // both negative
  if(sign_a && sign_b)
    return b.bits > a.bits;

  // a>=0 && b<=0
  return !sign_a && sign_b;
}
constexpr bool qge(const Float32 a, const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;

  const bool sign_a = a.bits >> Float32::SIGN_OFFSET;
  const bool sign_b = b.bits >> Float32::SIGN_OFFSET;

  // both positive
  if(!sign_a && !sign_b)
    return a.bits >= b.bits;
  // both negative
  if(sign_a && sign_b)
    return b.bits >= a.bits;

  // a>=0 && b<=0
  return !sign_a && sign_b;
}
constexpr bool qls(const Float32 a, const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;

  const bool sign_a = a.bits >> Float32::SIGN_OFFSET;
  const bool sign_b = b.bits >> Float32::SIGN_OFFSET;

  // both positive
  if(!sign_a && !sign_b)
    return a.bits < b.bits;
  // both negative
  if(sign_a && sign_b)
    return b.bits > a.bits;

  // a<=0 && b>=0
  return sign_a && !sign_b;
}
constexpr bool qle(const Float32 a, const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;

  const bool sign_a = a.bits >> Float32::SIGN_OFFSET;
  const bool sign_b = b.bits >> Float32::SIGN_OFFSET;

  // both positive
  if(!sign_a && !sign_b)
    return a.bits <= b.bits;
  // both negative
  if(sign_a && sign_b)
    return b.bits >= a.bits;

  // a<=0 && b>=0
  return sign_a && !sign_b;
}

/*
  5.10 
*/
constexpr bool total_order(const Float32 a, const Float32 b){
	//d)
	if(is_nan(a))
	{
		//3)
		if(is_nan(b))
		{
			//i)
			if(is_sign_minus(a))
				return true;

			//ii)
			//???

			//iii)
			//?????

			return false;
		}
		//1)
		else
		{
			//1)
			if(is_sign_minus(a))
				return true;
	
			return false;
		}
	}
	if(is_nan(b))
	{
		//2)
		if(is_sign_minus(b))
			return false;

		return true;
	}


	//a)
	if(qls(a, b))
		return true;
	//b)
	if(qgr(a, b))
		return false;

	//c)
	if(a == b)
		return true;
	//1) and 2)
	//if(is_Zero(a) && is_Zero(b)) ; not needed since if qeq() == true but a != b then they must be zeros of different sign
	return is_sign_minus(a);
}

constexpr bool total_order_mag(const Float32 a, const Float32 b){
  return total_order(abs(a), abs(b));
}


constexpr Float32 
    Float32::ZERO = Float32{0},
    Float32::NEG_ZERO = Float32{1u<<(Float32::SIGN_OFFSET)},
    Float32::INFTY{pack(IEEEFields{false, bit_mask<int32_t>(EXPONENT_SIZE), 0})},
    Float32::NEG_INFTY{neg(Float32::INFTY)},
    Float32::NaN{pack(IEEEFields{false, bit_mask<int32_t>(EXPONENT_SIZE), 1})},
    Float32::NEG_NaN{pack(IEEEFields{true, bit_mask<int32_t>(EXPONENT_SIZE), 1})},
    Float32::MAX{Float32::INFTY.bits - 1},
    Float32::MIN{Float32::NEG_INFTY.bits - 1};


int main() {
  // std::cout << "infinity test:\n" << std::fixed << std::setprecision(30)
  //           << "inf-1: " <<  std::bit_cast<float>(std::bit_cast<uint32_t>(std::numeric_limits<float>::infinity())-1) << '\n'
  //           << "max: " << std::numeric_limits<float>::max()  << '\n'
  //           << "equal: " << std::boolalpha << (std::bit_cast<float>(std::bit_cast<uint32_t>(std::numeric_limits<float>::infinity())-1) == std::numeric_limits<float>::max())
  //           << std::endl;

  float TEST_FLOATS[] = {0.1f, 0.5f, 1.75f, 10.10f, 0.001f, std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::max(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::min(), -0.0f, 0.0f};
  
  // test pack/unpack
  for(float f : TEST_FLOATS){
    Float32 my_float{std::bit_cast<uint32_t>(f)};
    Float32 unpacked_packed{pack(unpack(my_float))};

    assert(my_float==unpacked_packed && "pack/unpack identity");
  }

  // test binary operations
  for(auto a : TEST_FLOATS){
  for(auto b : TEST_FLOATS){
    Float32 af{static_cast<uint32_t>(a)};
    Float32 bf{static_cast<uint32_t>(b)};

    float a_mul_b = a*b;
    Float32 af_mul_bf = mul(af,bf);
    Float32 proper_mul{static_cast<uint32_t>(a_mul_b)};

    ASSERT_EQ(proper_mul, af_mul_bf, "Wrong result multiplication\n");
      // std::cout << "MyImpl: mantissa: " << unpack(af_mul_bf).mantissa << " exponent: " << unpack(af_mul_bf).exponent << std::endl;
      // std::cout << "BuiltIn: mantissa: " << unpack(proper_mul).mantissa << " exponent: " << unpack(proper_mul).exponent << std::endl;
  }}
}

