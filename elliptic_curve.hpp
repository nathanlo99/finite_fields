
#include "field.hpp"

#include <utility>

template <typename Field> struct EllipticCurve {
  using value_t = typename Field::value_t;
  using field_element_t = FieldElement<Field>;

  struct Point {
    field_element_t x, y;
    bool affine;
    Point(const Field *field)
        : x(field->value(0)), y(field->value(1)), affine(false) {}
    Point(const Field *field, const field_element_t x, const field_element_t y)
        : x(x), y(y), affine(true) {}
    static Point infinity(const Field *field) { return Point(field); }

    constexpr inline bool operator==(const Point &other) const {
      return x == other.x && y == other.y && affine == other.affine;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &point) {
      return os << "(" << point.x << " : " << point.y << " : " << point.affine
                << ")";
    }
  };

  using point_t = Point;

  const Field *field; // non-owning reference
  const field_element_t a, b;

  EllipticCurve(const Field *field, const field_element_t a,
                const field_element_t b)
      : field(field), a(a), b(b) {}
  EllipticCurve(const Field *field, const int a, const int b)
      : field(field), a(field->value(a)), b(field->value(b)) {}

  constexpr inline field_element_t f(const field_element_t x) const {
    return x * x * x + a * x + b;
  }
  constexpr inline Point affine_point(const field_element_t x,
                                      const field_element_t y) const {
    return Point(field, x, y);
  }
  constexpr inline Point affine_point(const int x, const int y) const {
    return affine_point(field->value(x), field->value(y));
  }

  constexpr inline bool is_point_on_curve(const point_t &point) const {
    if (!point.affine)
      return true;
    return point.y * point.y == f(point.x);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const EllipticCurve &curve) {
    return os << "Elliptic curve defined by y^2 = x^3 + " << curve.a.value
              << "*x + " << curve.b.value << " over " << *curve.field;
  }
};
