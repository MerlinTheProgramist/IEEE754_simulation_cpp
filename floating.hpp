#include <bitset>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <ostream>
#include <type_traits>

enum class RoundingType{
  ToEven, TowardZero, Positive, Negative, ToAway
};

struct Float32;
struct IEEEFields{
  bool sign_bit;
  // biased exponent
  uint16_t exponent;
  // mantissa with leading 1 (no subnormals)
  uint32_t mantissa;

  friend std::ostream& operator<<(std::ostream& os, const IEEEFields& f){
    os << "sign: " << f.sign_bit << " exponent: " << f.exponent << "(0b" << std::bitset<12>(f.exponent) << ")" << ", mantissa: " << f.mantissa << "(0b" << std::bitset<24>(f.mantissa) << ")";
    return os;
  }
};

template<typename T> 
constexpr T inline bit_mask(const size_t bits){
  return (1<<bits)-1;
}

/* 
  internal calculations
*/
constexpr IEEEFields unpack(Float32 f);
constexpr Float32 pack(IEEEFields fields);

/* 
  general operations
  section 5.3.1
*/

constexpr Float32 next_up(const Float32);
constexpr Float32 next_down(const Float32);
constexpr Float32 min_val(const Float32, const Float32);
constexpr Float32 max_val(const Float32, const Float32);
constexpr Float32 reminder(const Float32, const Float32);
/* 
  5.3.3 logBFormat operations
*/
constexpr Float32 scaleB(const Float32 x, const int32_t n);
constexpr int32_t logB(const Float32 x);

/*
  5.4.1 Arthmetic operations
*/
constexpr Float32 add(const Float32, const Float32);
constexpr Float32 sub(const Float32, const Float32);
constexpr Float32 mul(const Float32, const Float32);
constexpr Float32 div(const Float32, const Float32);
constexpr Float32 sqrt(const Float32);
constexpr Float32 fused_mul_add(const Float32, const Float32, const Float32);

/*
  5.4.2 conversion between float and string
*/
constexpr Float32 convert_from_decimal_character(std::string_view);
constexpr std::string convert_to_decimal_character(const Float32);
constexpr Float32 convert_from_hex_character(std::string_view);
constexpr std::string convert_to_hex_character(const Float32);

/*
  5.5.1 sign bit operations
*/
constexpr Float32 neg(const Float32);
constexpr Float32 abs(const Float32);
constexpr Float32 copy_sign(const Float32, const Float32);

/*
  5.6.1 comparisons
*/
constexpr bool qeq(const Float32, const Float32);
constexpr bool qne(const Float32, const Float32);
constexpr bool qgr(const Float32, const Float32);
constexpr bool qge(const Float32, const Float32);
constexpr bool qls(const Float32, const Float32);
constexpr bool qle(const Float32, const Float32);

/*
  5.7.1 non-computational operations
*/
constexpr bool is_754_version_1985 () {return false;}
constexpr bool is_754_version_2008 () {return false;}
constexpr bool is_754_version_2019 () {return true;}
constexpr uint32_t radix(){return 2;}

/* 
  5.7.2 general operations
*/
enum class cls {sNan, qNan, nInf, nSbn, nNrm, nZr, pZr, pSbn, pNrm, pInf};
constexpr cls  float_class(const Float32);
constexpr bool is_sign_minus(const Float32);
constexpr bool is_normal(const Float32);
constexpr bool is_finite(const Float32);
constexpr bool is_zero(const Float32);
constexpr bool is_subnormal(const Float32);
constexpr bool is_inf(const Float32);
constexpr bool is_nan(const Float32);
constexpr bool is_signaling(const Float32);
constexpr bool total_order(const Float32, const Float32);
constexpr bool total_order_mag(const Float32, const Float32);

struct Float32{
  uint32_t bits;

  static const uint32_t EXPONENT_SIZE=8, MANTISSA_SIZE=23;
  static const uint32_t EXPONENT_BIAS=127;
  static const uint32_t EXPONENT_FIRST_BIT=1<<(EXPONENT_SIZE+2);

  static const uint32_t TOTAL_SIZE = 1+EXPONENT_SIZE+MANTISSA_SIZE;
  static const uint32_t MANTISSA_OFFSET = 0;
  static const uint32_t EXPONENT_OFFSET = MANTISSA_SIZE;
  static const uint32_t SIGN_OFFSET = EXPONENT_SIZE + MANTISSA_SIZE;

  // implementation details
  static const uint32_t BIAS = 1<<(EXPONENT_SIZE+2); // 1024
  static const uint64_t MANTISSA_LEADING_1 = 1 << MANTISSA_SIZE; 
  static const uint32_t EXPONENT_MAX_BORDER = 1 << (13-1);

  static const Float32 ZERO, NEG_ZERO, INFTY, NEG_INFTY, QNaN, SNaN, MAX, MIN;

  constexpr bool operator==(const Float32& other)const{
    return this->bits==other.bits;
  }
  constexpr bool operator!=(Float32& other)const{
    return !(*this==other);
  }
  constexpr bool get_sign_bit() const{
    return (bool)((bits >> Float32::SIGN_OFFSET) & 1);
  }

  constexpr uint16_t get_exponent() const{
    return (static_cast<uint16_t>(bits >> Float32::EXPONENT_OFFSET) & bit_mask<uint16_t>(Float32::EXPONENT_SIZE))
           + Float32::BIAS - Float32::EXPONENT_BIAS;
  }
  constexpr uint32_t get_mantissa() const{
    return ((bits >> Float32::MANTISSA_OFFSET) & bit_mask<uint32_t>(Float32::MANTISSA_SIZE));
  }
  friend std::ostream& operator<<(std::ostream& os, const Float32& f){
    os << std::bit_cast<float>(f.bits) << " (0b" << std::bitset<TOTAL_SIZE>(f.bits) << ")";
    return os;
  }
  constexpr operator float() const{
    return std::bit_cast<float>(this->bits);
  }
  constexpr Float32(const uint32_t bits)
    :bits{bits}
  {}
};



constexpr IEEEFields unpack(Float32 f){
  auto fields = IEEEFields{
    .sign_bit=f.get_sign_bit(),
    .exponent=f.get_exponent(),
    .mantissa=f.get_mantissa()
  };

  // if exponent equals zero
  if(fields.exponent == Float32::BIAS - Float32::EXPONENT_BIAS){
    if(fields.mantissa == 0){
      fields.exponent = 0;
      return fields;
    }

    fields.mantissa <<= 1;
    // change subnormals to normals
    while(fields.mantissa < Float32::MANTISSA_LEADING_1){
      fields.mantissa <<= 1;
      fields.exponent -= 1;
    }
  }

  fields.mantissa |= Float32::MANTISSA_LEADING_1;
  return fields;
}

constexpr Float32 pack(const IEEEFields fields){
  	uint32_t sign_bit = fields.sign_bit;
  	uint32_t exponent = fields.exponent;
  	uint32_t mantissa = fields.mantissa;

  	//arbitrary border for overflow detection
  	if(exponent > Float32::EXPONENT_MAX_BORDER)
  		return Float32{(sign_bit << Float32::SIGN_OFFSET) | 0};
  	
  	if(fields.exponent == (bit_mask<uint16_t>(Float32::EXPONENT_SIZE) + static_cast<uint16_t>(Float32::BIAS) - static_cast<uint16_t>(Float32::EXPONENT_BIAS)))
  	  return ((mantissa&bit_mask<uint32_t>(Float32::MANTISSA_SIZE))!=0) 
              ? Float32::QNaN
  	          : ((sign_bit)
  	            ? Float32::NEG_INFTY
  	            : Float32::INFTY);

  	if(exponent > Float32::BIAS + Float32::EXPONENT_BIAS)
  		return sign_bit ? Float32::MIN : Float32::MAX;

  	if(exponent <= Float32::BIAS - Float32::EXPONENT_BIAS)
  		exponent--;

  	//convert from extended precision normal to subnormal
  	while(exponent < Float32::BIAS - Float32::EXPONENT_BIAS
  	   && mantissa != 0)
  	{
  		exponent++;
  		mantissa >>= 1;
  	}

  	//if no implicit one, return 0 because extended precision 
  	if(mantissa == 0) 
  		return {sign_bit << Float32::SIGN_OFFSET};
	
  	return {((sign_bit & 1)  << Float32::SIGN_OFFSET)
  	     | (((exponent - Float32::BIAS + Float32::EXPONENT_BIAS) & bit_mask<uint32_t>(Float32::EXPONENT_SIZE)) << Float32::EXPONENT_OFFSET)
  	     | (  mantissa & bit_mask<uint32_t>(Float32::MANTISSA_SIZE))};
}

static_assert(Float32::TOTAL_SIZE == sizeof(Float32)*8);


template<typename T, const RoundingType=RoundingType::TowardZero>
constexpr Float32 convert_from_int(T i){
  static_assert(std::is_integral<T>() && std::is_signed<T>() && "Converion type must be signed integral");
  using uT = std::make_unsigned_t<T>;

  if(i == 0) 
    return Float32::ZERO;
  
  bool sign = i<0;
  i = std::abs(i);
  
  uT mantissa = static_cast<uT>(i);
  uint32_t exponent = Float32::EXPONENT_BIAS + Float32::MANTISSA_SIZE;

  static_assert(static_cast<T>(Float32::MANTISSA_LEADING_1) < bit_mask<T>(Float32::MANTISSA_SIZE+1));
  while(mantissa > bit_mask<T>(Float32::MANTISSA_SIZE+1)){
    mantissa >>= 1;
    exponent += 1;
  }
    while(mantissa < static_cast<T>(Float32::MANTISSA_LEADING_1)){
      mantissa <<= 1;
      exponent -= 1;
  }

  exponent <<= Float32::EXPONENT_OFFSET;
  mantissa &= bit_mask<T>(Float32::MANTISSA_SIZE);

  return Float32{sign<<Float32::SIGN_OFFSET | static_cast<uint32_t>(exponent) | static_cast<uint32_t>(mantissa)};
}

template<typename T, const RoundingType=RoundingType::TowardZero>
constexpr T round_to_int(const Float32 fa){
  static_assert(std::is_integral<T>() && std::is_signed<T>() && "Converion type must be signed integral");

  constexpr size_t INT_SIZE = 8*sizeof(T);
  if(is_nan(fa))
    return std::numeric_limits<T>::min();

  if(is_inf(fa)){
    return is_sign_minus(fa)
         ? std::numeric_limits<T>::min()
         : std::numeric_limits<T>::min();
  }

  auto a = unpack(fa);
  T mantissa = static_cast<T>(a.mantissa);

  if(a.exponent < Float32::BIAS)
    mantissa = 0;
  else if(a.exponent >= Float32::BIAS + INT_SIZE-1)
    return a.sign_bit
         ? std::numeric_limits<T>::min()
         : std::numeric_limits<T>::max();

  if(a.exponent > Float32::BIAS + Float32::MANTISSA_SIZE)
    mantissa <<= a.exponent - Float32::BIAS - Float32::MANTISSA_SIZE;
  else
    mantissa >>= Float32::MANTISSA_SIZE - (a.exponent - Float32::BIAS);

  return a.sign_bit ? -(mantissa) : mantissa;
}

