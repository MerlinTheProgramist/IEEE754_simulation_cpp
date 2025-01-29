#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

#include "floating.hpp"

template<typename T> 
constexpr bool EQ (const T a, const T b){return a==b;}
constexpr bool EQ (const float a, const float b){
  return std::bit_cast<uint32_t>(a)==std::bit_cast<uint32_t>(b) 
      || ((std::abs(a - b) <= std::numeric_limits<float>::epsilon() * std::max(std::abs(a), std::abs(b))))
      || ((std::abs(b - a) <= std::numeric_limits<float>::epsilon() * std::max(std::abs(a), std::abs(b))));
}
constexpr bool EQ (const Float32 a, const Float32 b){
  return EQ((float)a, (float)b);
}

#define RETURN_IF_A_B_BOTH(a, b, only_a, only_b, both) \
  if((a)){ \
    if((b)) \
      return both; \
    return (only_a); \
  }else if((b)){ \
    return (only_b); \
  }

constexpr bool is_nan(const Float32 a){
  // max exponent && mantissa != 0
  return ((a.bits >> Float32::EXPONENT_OFFSET) & bit_mask<uint32_t>(Float32::EXPONENT_SIZE)) == bit_mask<uint32_t>(Float32::EXPONENT_SIZE)
      && (a.bits & bit_mask<uint32_t>(Float32::MANTISSA_SIZE)) != 0;
}
constexpr bool is_signaling(const Float32 a){
  return a == Float32::SNaN;
}
constexpr bool is_inf(Float32 a){
  // max exponent && mantissa = 0
  return ((a.bits >> Float32::EXPONENT_OFFSET) & bit_mask<uint32_t>(Float32::EXPONENT_SIZE)) == bit_mask<uint32_t>(Float32::EXPONENT_SIZE)
      && (a.bits & (bit_mask<uint32_t>(Float32::MANTISSA_SIZE))) == 0;
}
constexpr bool is_zero(Float32 a){
  return (a.bits & bit_mask<uint32_t>(Float32::MANTISSA_SIZE + Float32::EXPONENT_SIZE)) == 0;
}
constexpr bool is_sign_minus(const Float32 a){
  return a.bits&(1<<Float32::SIGN_OFFSET);
}
constexpr bool is_normal(const Float32 a){
  if(is_nan(a) || is_inf(a))
    return false;
  return !is_subnormal(a);
}
constexpr bool is_subnormal(const Float32 a){
  return ((a.bits >> Float32::EXPONENT_OFFSET) & bit_mask<uint32_t>(Float32::EXPONENT_SIZE)) == 0;
}
constexpr bool is_finite(const Float32 a){
  return !is_inf(a);
}
constexpr bool total_order_mag(const Float32 a, const Float32 b){
  return total_order(abs(a), abs(b));
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

    auto a = unpack(fa);
    auto b = unpack(fb);

    const uint8_t sign = 0;

	//adjust smaller number and simply add their mantisa
	//then adjust result
	uint16_t exponent{};
	uint32_t to_shift{};
	uint32_t mantissa{};
	if(a.exponent < b.exponent)
	{
		exponent = b.exponent; 
		to_shift = b.exponent - a.exponent;
		mantissa  = (a.mantissa >> to_shift) + b.mantissa;
		if(to_shift >= Float32::TOTAL_SIZE)
			mantissa = b.mantissa;
	}
	else
	{
		exponent = a.exponent; 
		to_shift = a.exponent - b.exponent;
		mantissa  = (b.mantissa >> to_shift) + a.mantissa;
		if(to_shift >= 32)
			mantissa = a.mantissa;
	}
	while(mantissa > bit_mask<uint32_t>(Float32::MANTISSA_SIZE+1)) //24 bits
	{
		exponent++;
		mantissa >>= 1;
	} 

	const Float32 t = pack({sign, exponent, mantissa});
	if(is_nan(t))
		return sign ? Float32::MIN : Float32::MAX;

	return t;
    
}


constexpr Float32 sub(const Float32 fa, const Float32 fb){
	if(is_nan(fa))
		return fa;
	if(is_nan(fb))
		return fb;

	if(is_sign_minus(fb))
		return add(fa, neg(fb));
	if(is_sign_minus(fa))
		return neg(add(neg(fa), fb));


	if(is_zero(fa)){
	  if(is_zero(fb))
	    return Float32::ZERO;
	  return neg(fb);
	}
	if(is_zero(fb))
	  return fa;

	if(is_inf(fa)){
	  if(is_inf(fb))
	    return neg(Float32::QNaN);
	  return Float32::INFTY;
	}
	if(is_inf(fb))
	  return Float32::NEG_INFTY;
	

	auto a = unpack(fa);
	auto b = unpack(fb);

	const bool sign = qls(fa, fb);

	uint16_t exponent{}; 
	uint32_t to_shift{}; 
	uint64_t a_mantisa = a.mantissa; 
	uint64_t b_mantisa = b.mantissa; 
	uint64_t mantissa{};

  const size_t PREC_SHIFT = 40;
	
	a_mantisa <<= PREC_SHIFT;
	b_mantisa <<= PREC_SHIFT;

	if(!sign)
	{
		exponent = a.exponent;
		to_shift = a.exponent - b.exponent;
		mantissa  = a_mantisa - (b_mantisa >> to_shift);

		if(to_shift >= sizeof(mantissa)*8)
			mantissa = a_mantisa - 1;
		
	}
	else
	{
		exponent = b.exponent;
		to_shift = b.exponent - a.exponent;
		mantissa  = b_mantisa - (a_mantisa >> to_shift);
		if(to_shift >= sizeof(mantissa)*8)
			mantissa = b_mantisa - 1;
	}

	if(mantissa == 0)
		return 0;

	while(mantissa > (bit_mask<uint64_t>(Float32::MANTISSA_SIZE+1) << PREC_SHIFT)) //24 bits
	{
		exponent++;
		mantissa >>= 1;
	} 
	
	while(mantissa < (Float32::MANTISSA_LEADING_1 << PREC_SHIFT)) //24 bits
	{
		exponent--;
		mantissa <<= 1;
	} 
	
	if(mantissa == 0)
		return 0;
	mantissa >>= PREC_SHIFT;
	
	const Float32 t = pack({sign,  exponent, static_cast<uint32_t>(mantissa)});
	if(is_nan(t))
		return sign ? Float32::MIN : Float32::MAX;

	return t;
} 

constexpr Float32 div(const Float32 fa, const Float32 fb){
	if(is_nan(fa))
		return fa;
	if(is_nan(fb))
		return fb;
	
	const bool sign = is_sign_minus(fa) ^ is_sign_minus(fb);
	
	RETURN_IF_A_B_BOTH(is_inf(fa), is_inf(fb), 
      	            sign?Float32::NEG_INFTY:Float32::INFTY, 
      	            sign?Float32::NEG_ZERO:Float32::ZERO, 
      	            neg(Float32::QNaN))
	RETURN_IF_A_B_BOTH(is_zero(fa), is_zero(fb), 
      	            sign?Float32::NEG_ZERO:Float32::ZERO, 
      	            sign?Float32::NEG_INFTY:Float32::INFTY,
      	            neg(Float32::QNaN))
	
	IEEEFields a = unpack(fa);
	IEEEFields b = unpack(fb);

	a.sign_bit ^= b.sign_bit;

	/*
		shifted 40 to increase precision of result
		later res is shifted 40 to the right
	*/
  const size_t PREC_SHIFT = 40;
	uint64_t num = static_cast<uint64_t>(a.mantissa) << PREC_SHIFT;
	uint64_t div = b.mantissa;
	uint64_t res = num / div;

	res <<= Float32::MANTISSA_SIZE;

	/*
		subtracting exponents so BIAS has to be added
	*/
	a.exponent = a.exponent - b.exponent + Float32::BIAS;

	if(res == 0)
		return a.sign_bit ? Float32::MIN : Float32::MAX;

	while((res & (1ul<<(8*sizeof(res)-1))) == 0)
	{
		res <<= 1;
		a.exponent--;
	}
	
	a.mantissa = res >> PREC_SHIFT;
	const Float32 result = pack(a);
	if(is_nan(result))
		return a.sign_bit?Float32::MIN:Float32::MAX;

	return result;
  
}
constexpr Float32 sqrt(const Float32 fa){
  if(is_zero(fa) || is_nan(fa))
    return fa;
  if(is_sign_minus(fa))
    return neg(Float32::QNaN);
  if(is_inf(fa)) // negative 0 is nan
    return fa;

  auto a = unpack(fa);

  /*
	sqrt(mant) can be computed as integer

	in case of odd number we have to multiply mantisa by 2 first, then take square root
	because 2^(odd / 2) will lead to 2^(even + 1/2) = 2^even * sqrt(2) 
	sqrt(2 * mantisa)  = sqrt(mantisa) * sqrt(2)
	*/

	a.exponent -= Float32::BIAS;
	const bool even = a.exponent%2==0;
	a.exponent = static_cast<int16_t>(a.exponent) >> 1;

	a.exponent += Float32::BIAS;

	uint64_t mantissa = a.mantissa;
	mantissa = std::sqrt(mantissa << (40 - even));
	mantissa >>= Float32::EXPONENT_SIZE;
	a.mantissa = mantissa;

	return pack(a);
}

constexpr Float32 fused_mul_add(const Float32 x, const Float32 y, const Float32 z){
  return add(mul(x,y), z);
}

constexpr cls float_class(const Float32 a){
  using enum cls;
  if(is_nan(a))
    return qNan;
  if(is_signaling(a))
    return sNan;
  if(is_inf(a))
    return is_sign_minus(a)?nInf:pInf;
  if(is_zero(a))
    return is_sign_minus(a)?nZr:pZr;
  if(is_normal(a))
    return is_sign_minus(a)?nNrm:pNrm;
  return is_sign_minus(a)?nSbn:pSbn;
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
      return neg(Float32::QNaN);

    if(is_zero(fa))
      return fa;
    if(is_zero(fb))
      return fb;
    
    // sign remains only if sign is diffrent
    const bool sign_bit = (fa.bits ^ fb.bits) >> Float32::SIGN_OFFSET;
    if(is_inf(fa) || is_inf(fb))
      return Float32{(static_cast<uint32_t>(sign_bit) << Float32::SIGN_OFFSET) | Float32::INFTY.bits};
    
    IEEEFields a = unpack(fa);
    IEEEFields b = unpack(fb);
        
    /*
      exponent are added
      mantissas are mulilied and lower bits are discarded via shift
    */
      
    uint16_t exponent = a.exponent + b.exponent - Float32::BIAS;
    uint64_t mantissa = static_cast<uint64_t>(a.mantissa)
                      * static_cast<uint64_t>(b.mantissa);

    mantissa >>= Float32::MANTISSA_SIZE;

    // normalization
    if(mantissa != 0){
      // check if mantissa is greater than its maximum
      while(mantissa > bit_mask<uint64_t>(Float32::MANTISSA_SIZE+1)){
        exponent += 1;
        mantissa >>= 1;
      }
      // check if mantissa is less than its minimum
      while(mantissa < Float32::MANTISSA_LEADING_1){
        exponent -= 1;
        mantissa <<= 1;
      }
    }
    mantissa += 1;

    const Float32 result = pack(IEEEFields{sign_bit,exponent,static_cast<uint32_t>(mantissa)});

    // prevent nan
    if(is_nan(result))
      return sign_bit ? Float32::MIN : Float32::MAX;

    return result;
}



constexpr Float32 scaleB(const Float32 x, const int32_t N){
  if(N==0 || is_inf(x) || is_nan(x) || is_zero(x))
    return x;

  if(N >= static_cast<int32_t>(Float32::BIAS))
    return Float32::MAX;
  if(N <= -static_cast<int32_t>(Float32::BIAS))
    return Float32::MIN;
  
  auto a = unpack(x);
  a.exponent += N;
  auto f = pack(a);
  return is_inf(f) ? Float32::MAX : f;
}

constexpr Float32 neg(const Float32 x){
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
  if(is_nan(a) || is_nan(b)
  ||(is_zero(a) && is_zero(b)))
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
  if(is_zero(a) && is_zero(b))
    return true;

  const bool sign_a = is_sign_minus(a);
  const bool sign_b = is_sign_minus(b);

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
  if(is_nan(a) || is_nan(b) 
  ||(is_zero(a) && is_zero(b)))
    return false;

  const bool sign_a = a.bits >> Float32::SIGN_OFFSET;
  const bool sign_b = b.bits >> Float32::SIGN_OFFSET;

  // both positive
  if(!sign_a && !sign_b)
    return a.bits < b.bits;
  // both negative
  if(sign_a && sign_b)
    return b.bits < a.bits;

  // a<=0 && b>=0
  return sign_a && !sign_b;
}
constexpr bool qle(const Float32 a, const Float32 b){
  if(is_nan(a) || is_nan(b))
    return false;
  if(is_zero(a) && is_zero(b))
    return true;

  const bool sign_a = a.bits >> Float32::SIGN_OFFSET;
  const bool sign_b = b.bits >> Float32::SIGN_OFFSET;

  // both positive
  if(!sign_a && !sign_b)
    return a.bits <= b.bits;
  // both negative
  if(sign_a && sign_b)
    return b.bits <= a.bits;

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

constexpr Float32 next_up(const Float32 a){
  if(is_nan(a))
    return a;
  if(is_zero(a))
    return 1;
  if(is_sign_minus(a))
    return {a.bits-1};
  else 
    return is_inf(a) ? a : Float32{a.bits+1};
}
constexpr Float32 next_down(const Float32 a){
  return neg(next_up(neg(a)));
}
constexpr Float32 min_val(const Float32 a, const Float32 b){
  if(is_nan(a))
    return a;
  if(is_nan(b))
    return b;

  return qls(a,b)? a : b;
}
constexpr Float32 max_val(const Float32 a, const Float32 b){
  if(is_nan(a))
    return a;
  if(is_nan(b))
    return b;

  return qgr(a,b)? a : b;
}
// constexpr Float32 reminder(const Float32, const Float32){}



constexpr Float32 
    Float32::ZERO = Float32{0},
    Float32::NEG_ZERO = Float32{1u<<Float32::SIGN_OFFSET},
    Float32::INFTY{bit_mask<int32_t>(EXPONENT_SIZE) << (Float32::EXPONENT_OFFSET)},
    Float32::NEG_INFTY{neg(Float32::INFTY)},
    Float32::QNaN{(bit_mask<uint32_t>(EXPONENT_SIZE) << (Float32::EXPONENT_OFFSET)) | (1<<(MANTISSA_SIZE-1))},
    Float32::SNaN{(bit_mask<uint32_t>(EXPONENT_SIZE) << (Float32::EXPONENT_OFFSET)) | (1<<(MANTISSA_SIZE-2))},
    Float32::MAX{Float32::INFTY.bits - 1},
    Float32::MIN{neg(Float32::MAX)};


constexpr bool operator==(const float a, const Float32 b){
  return std::bit_cast<uint32_t>(a)== b.bits;
}
constexpr bool operator!=(const float a, const Float32 b){
  return !(a == b);
}




