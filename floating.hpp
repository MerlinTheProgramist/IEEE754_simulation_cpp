#include <bitset>
#include <cstddef>
#include <cstdint>
#include <ostream>

enum class RoundingType{
  ToEven, TowardZero, Positive, Negative, ToAway
};

struct Float32;
struct IEEEFields{
  bool sign_bit;
  // biased exponent
  int16_t exponent;
  // mantissa with leading 1 (no subnormals)
  uint32_t mantissa;
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

template<typename T, const RoundingType>
constexpr T round_to_integral(const Float32);
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
template<typename T, const RoundingType>
constexpr Float32 convert_from_int(T);

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

/* 
  5.7.2 general operations
*/
enum class cls {sNan, qNan, nInf, nNrm, nZr, pZr, pSbn, pNrm, pInf};
constexpr bool is_sign_minus(const Float32);
constexpr bool is_normal(const Float32);
constexpr bool is_finite(const Float32);
constexpr bool is_zero(const Float32);
constexpr bool is_subnormal(const Float32);
constexpr bool is_infinite(const Float32);
constexpr bool is_nan(const Float32);
constexpr bool is_signaling(const Float32);
constexpr bool is_canonical(const Float32);
constexpr bool total_order(const Float32, const Float32);
constexpr bool total_order_mag(const Float32, const Float32);

struct Float32{
  uint32_t bits;

  static const uint32_t BIAS = 1024;
  static const uint32_t EXPONENT_SIZE=8, MANTISSA_SIZE=23;
  static const int32_t EXPONENT_BIAS=127;
  static const int32_t EXPONENT_FIRST_BIT=1<<(EXPONENT_SIZE+2);

  static const uint32_t TOTAL_SIZE = 1+EXPONENT_SIZE+MANTISSA_SIZE;
  static const uint32_t MANTISSA_OFFSET = 0;
  static const uint32_t EXPONENT_OFFSET = MANTISSA_SIZE;
  static const uint32_t SIGN_OFFSET = EXPONENT_SIZE + MANTISSA_SIZE;

  static const Float32 ZERO, NEG_ZERO, INFTY, NEG_INFTY, NaN, NEG_NaN, MAX, MIN;

  constexpr bool operator==(const Float32& other)const{
    return this->bits==other.bits;
  }
  constexpr bool operator!=(Float32& other)const{
    return !(*this==other);
  }
  constexpr bool get_sign_bit() const{
    return (bool)((bits >> Float32::SIGN_OFFSET) & 1);
  }

  constexpr int16_t get_exponent() const{
    return static_cast<int16_t>(
              (bits >> Float32::EXPONENT_OFFSET) & bit_mask<uint16_t>(Float32::EXPONENT_SIZE)
            ) - Float32::EXPONENT_BIAS;
  }
  constexpr uint32_t get_mantissa() const{
    return ((bits >> Float32::MANTISSA_OFFSET) & bit_mask<uint32_t>(Float32::MANTISSA_SIZE));
  }

  friend std::ostream& operator<<(std::ostream& os, const Float32& f){
    auto a = unpack(f);
    os << "bits: " << std::bitset<TOTAL_SIZE>(f.bits) << '\n'
       << "sign: " << a.sign_bit << "exponent: " << a.exponent << ", mantissa: " << a.mantissa;
    return os;
  }
  
};



constexpr IEEEFields unpack(Float32 f){
  return IEEEFields{
    .sign_bit=f.get_sign_bit(),
    .exponent=f.get_exponent(),
    .mantissa=f.get_mantissa()
  };
}

constexpr Float32 pack(IEEEFields fields){
  return Float32{.bits=(fields.sign_bit << (Float32::SIGN_OFFSET))
    | ((static_cast<uint32_t>(fields.exponent + Float32::EXPONENT_BIAS) & bit_mask<uint32_t>(Float32::EXPONENT_SIZE)) << Float32::EXPONENT_OFFSET)
    | ((fields.mantissa & bit_mask<uint32_t>(Float32::MANTISSA_SIZE)) << Float32::MANTISSA_OFFSET)
    };
}

static_assert(Float32::TOTAL_SIZE == sizeof(Float32)*8);




constexpr Float32 int_to_float(int32_t integer);
